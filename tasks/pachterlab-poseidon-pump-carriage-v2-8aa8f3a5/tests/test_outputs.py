"""Geometric checks on cylindrical seats and straight probe insertion."""

import ctypes
import hashlib
import json
import math
import os
import stat
import struct
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from functools import cache, cached_property
from pathlib import Path
from select import select
import shutil
import subprocess
import sys
import tempfile
from typing import NamedTuple

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
import shapely
from scipy.optimize import linear_sum_assignment
from shapely import MultiLineString

GRADER = Path(__file__).resolve().parent / "grader"
sys.path.insert(0, str(GRADER))
from clearance import access_rows
from lands import get_lands, wall_rows
from poses import FACET, repose


OUT = Path("/app/out")
INPUTS = Path("/app/inputs")
PROGRAM = OUT / "bore_access.py"
PYTHON = "/opt/conda/bin/python"
UNPRIVILEGED = dict(user=65534, group=65534, extra_groups=[])
SAVED = {
    "base": "Scope Base v33.stl",
    "screen-leg": "Scope Screen Leg Mount v33.stl",
    "top-stage": "Scope Top Stage v33.stl",
}
FRESH = ("pump-base", "carriage", "brace", "slide")
LISTED_CLI = {f"{key}-cli": key for key in SAVED}
POSED_SCREEN = ("screen-familiar-name", "screen-unrelated-name")
REPLAYS = ("screen-leg", "carriage", "slide")
FINAL = {f"{key}-after": key for key in SAVED}
ALL = (*SAVED, *LISTED_CLI, "screen-independent", *POSED_SCREEN, *FRESH,
       *(f"{key}-final-cli" for key in REPLAYS), *FINAL)
WITH_LANDS = tuple(key for key in ALL if key not in {"slide", "slide-final-cli"})
ALL_CASES = pytest.mark.parametrize("key", ALL)
LAND_CASES = pytest.mark.parametrize("key", WITH_LANDS)


def scratch_dir():
    path = Path(tempfile.mkdtemp(dir="/tmp"))
    os.chown(path, 65534, 65534)
    return path


def run_delivered(mesh, destination, scratch, *, cwd=None, model_only=False):
    assert PROGRAM.is_file(), f"Missing required program {PROGRAM}"
    command = [PYTHON, str(PROGRAM), "--mesh", str(mesh), "--out", str(destination)]
    env = {**os.environ, "HOME": str(scratch), "TMPDIR": str(scratch),
           "MPLCONFIGDIR": str(scratch)}
    env.pop("PYTEST_CURRENT_TEST", None)
    with tempfile.TemporaryFile() as out, tempfile.TemporaryFile() as err:
        options = ["--model-only"] if model_only else []
        return subprocess.run(
            [sys.executable, "-I", "-S", GRADER / "launcher.py", cwd or scratch, *options, *command],
            cwd="/", env=env, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
        ).returncode


# Inspect the file descriptor under the publisher's identity. Opening the
# resolved descriptor afterwards keeps the inode fixed while its bytes are read.
REPORT_READER = """
import os
from pathlib import Path
import shutil
import stat
import sys

path = Path(sys.argv[1])
root = Path(sys.argv[2])
directory = path.parent.resolve(strict=True)
if not directory.is_relative_to(root):
    sys.exit(1)
fd = os.open(path, os.O_PATH)
info = os.fstat(fd)
descriptor = Path(f"/proc/self/fd/{fd}")
target = descriptor.readlink()
if not stat.S_ISREG(info.st_mode) or not target.is_relative_to(root):
    sys.exit(1)
with descriptor.open("rb") as source:
    out = sys.stdout.buffer
    target_hex = os.fsencode(target).hex()
    out.write(f"{info.st_dev} {info.st_ino} {target_hex}\\n".encode("ascii"))
    shutil.copyfileobj(source, out)
"""


class FileSnapshot(NamedTuple):
    identity: tuple[int, int] | None
    contents: bytes | None
    target: Path | None


def regular_snapshot(path, root=None):
    """Read a regular file, its identity and resolved path under the publisher's identity."""
    if root is None:
        root = OUT if path.is_relative_to(OUT) else path.parent
    read = subprocess.run(
        [PYTHON, "-I", "-S", "-c", REPORT_READER, str(path), str(root)],
        cwd="/", stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, **UNPRIVILEGED,
    )
    if read.returncode:
        return FileSnapshot(None, None, None)
    header, content = read.stdout.split(b"\n", 1)
    dev, inode, target = header.decode("ascii").split()
    return FileSnapshot((int(dev), int(inode)), content, Path(os.fsdecode(bytes.fromhex(target))))


def regular_bytes(path, root=None):
    return regular_snapshot(path, root).contents


def read_triangles(path):
    with path.open("rb") as stream:
        stream.seek(84)
        return np.fromfile(stream, dtype=FACET)["vertices"].astype(float)


def unique_object(pairs):
    result = {}
    for name, value in pairs:
        if name in result:
            raise ValueError(f"Repeated JSON member {name!r}")
        result[name] = value
    return result


@dataclass
class Survey:
    key: str
    raw: bytes | None
    source: np.ndarray
    rotation: np.ndarray
    translation: np.ndarray
    triangles: np.ndarray
    status: int = 0

    @cached_property
    def document(self):
        assert self.status == 0, f"{self.key}: valid whole-component input exited {self.status}"
        assert self.raw is not None, f"{self.key}: survey.json missing from requested directory"
        try:
            result = json.loads(self.raw, object_pairs_hook=unique_object)
        except (ValueError, UnicodeError, RecursionError) as error:
            assert False, f"{self.key}: unreadable JSON ({error})"
        assert isinstance(result, dict), f"{self.key}: a survey must be a JSON object"
        return result

    @cached_property
    def bores(self):
        records = self.document.get("bores")
        assert isinstance(records, list), f"{self.key}: bores must be a JSON array"
        return records

    @cached_property
    def lands(self):
        placed_in_source_order = (
            (self.source @ self.rotation.T + self.translation).astype("<f4").astype(float)
        )
        return get_lands(self.key, self.source, self.rotation, self.translation,
                         placed_triangles=placed_in_source_order)

    @cached_property
    def pairs(self):
        records, lands = self.bores, self.lands
        assert len(records) == len(lands), (
            f"{self.key}: expected {len(lands)} qualifying cylinders, found {len(records)}"
        )
        c = vectors(records, "c")
        d = vectors(records, "d")
        radii = np.array([radius(b) for b in records])
        distances = np.linalg.norm(c[:, None] - [land.c for land in lands], axis=2)
        radial_errors = np.abs(radii[:, None] - [land.r for land in lands])
        angles = unsigned_angles(d, np.asarray([land.d for land in lands]))
        # A feasible perfect matching always wins over any matching with a
        # rejected pair. Within that choice, prefer the closest nominal axes.
        worst = np.maximum.reduce((distances / 0.03, radial_errors / 0.03, angles / 0.03))
        cost = 1.0 - 1.0 / (1.0 + worst) + (len(records) + 1) * (worst > 1.0)
        ii, jj = linear_sum_assignment(cost)
        return [(records[i], lands[j], distances[i, j], radial_errors[i, j], angles[i, j])
                for i, j in zip(ii, jj, strict=True)]

    @cached_property
    def scores(self):
        scores = []
        for report, land, *_ in self.pairs:
            rows = report.get("rows")
            check_rows(rows)
            d = land.d if np.dot(land.d, report["d"]) >= 0 else -land.d
            e = vector(report, "e")
            e -= np.dot(e, d) * d
            norm = np.linalg.norm(e)
            assert norm > 0, "Clock zero must project to a direction in the cylinder plane"
            e /= norm
            angles = np.deg2rad(np.arange(0, 360, 2))[:, None]
            normals = np.cos(angles) * e + np.sin(angles) * np.cross(d, e)
            walls = wall_rows(land, d, normals)
            reference = access_rows(self.triangles, land.c, d, e, land.r, walls)
            weight = sum(b - a for wall in walls for a, b in wall)
            assert weight > 0, f"No clock stations intersect reference cylinder at {land.c}"
            scores.append([column_difference(reference, rows, j) / weight for j in range(3)])
        return np.asarray(scores)


def is_number(value):
    return type(value) in (float, int) and math.isfinite(value)


def vector(record, key):
    assert isinstance(record, dict), "A bore record must be a JSON object"
    values = record.get(key)
    assert isinstance(values, list) and len(values) == 3 and all(map(is_number, values)), (
        f"{key} must be a finite XYZ triple"
    )
    return np.array(values, dtype=float)


def vectors(records, key):
    return np.array([vector(record, key) for record in records])


def radius(record):
    assert isinstance(record, dict) and is_number(record.get("r")), "r must be a finite radius"
    return record["r"]


def unsigned_angles(axes, nominal):
    cross = np.cross(axes[:, None, :], nominal[None, :, :])
    return np.rad2deg(np.arctan2(np.linalg.norm(cross, axis=2), np.abs(axes @ nominal.T)))


def check_rows(rows):
    assert isinstance(rows, list) and len(rows) == 180, "Expected rows for k = 0 through 179"
    for i, row in enumerate(rows):
        assert isinstance(row, list) and len(row) == 3, f"Row {i} needs wall, plus and minus"
        for entry in row:
            assert isinstance(entry, list), f"Row {i} columns must be interval lists"
            end = -math.inf
            for pair in entry:
                assert isinstance(pair, list) and len(pair) == 2 and all(map(is_number, pair)), (
                    f"Row {i} contains an invalid pair of stations"
                )
                a, b = pair
                assert end < a < b, f"Row {i} intervals must be positive, sorted and merged"
                end = b


def check_access_subsets(rows):
    check_rows(rows)
    roundoff = 1e-9  # mm, for comparisons of axial stations
    for i, (wall, plus, minus) in enumerate(rows):
        spans = []
        for low, high in wall:
            if spans and low - spans[-1][1] <= roundoff:
                spans[-1][1] = high
            else:
                spans.append([low, high])
        for name, access in (("plus", plus), ("minus", minus)):
            for start, end in access:
                assert any(low - roundoff <= start and end <= high + roundoff
                           for low, high in spans), (
                    f"Row {i} {name} interval [{start}, {end}] extends outside its wall"
                )


def column_difference(reference, reported, column):
    def lines(row):
        return MultiLineString([[(a, 0), (b, 0)] for a, b in row[column]])

    ref = [lines(row) for row in reference]
    got = [lines(row) for row in reported]
    return shapely.length(shapely.symmetric_difference(ref, got)).sum()


@cache
def geometry(component, posed):
    path = INPUTS / SAVED[component] if component in SAVED else GRADER / "meshes" / f"{component}.stl"
    original = read_triangles(path)
    if posed:
        with tempfile.TemporaryDirectory() as work:
            mesh = Path(work) / "component.stl"
            triangles, rotation, translation = repose(component, path, mesh)
            content = mesh.read_bytes()
    else:
        triangles, rotation, translation = original, np.eye(3), np.zeros(3)
        content = path.read_bytes()
    return content, (original, rotation, translation, triangles)


class Inspections:
    """Keep each distinct report with the mesh and placement it must describe."""

    def __init__(self, inventory):
        self.inventory = inventory
        self.groups = {key: {} for key in ALL}
        self.observations = {}
        self.failures = {}
        for key in SAVED:
            self.record(key, key, False, regular_bytes(OUT / key / "survey.json"))

    def record(self, group, component, posed, raw, status=0):
        identity = (component, posed, raw, status)
        if identity not in self.observations:
            _, placement = geometry(component, posed)
            self.observations[identity] = Survey(component, raw, *placement, status)
        self.groups[group][identity] = self.observations[identity]

    def saved_checkpoint(self):
        for component in SAVED:
            self.record(f"{component}-after", component, False,
                        regular_bytes(OUT / component / "survey.json"))

    def invoke(self, group, component, *, posed=False, name=None, destination=None,
               model_only=False, published=False, relative=False, listed_path=False):
        scratch = scratch_dir()
        try:
            caller = scratch / "caller"
            caller.mkdir()
            os.chown(caller, 65534, 65534)
            if published or listed_path:
                mesh = INPUTS / SAVED[component]
            else:
                source_dir = scratch / "parts"
                source_dir.mkdir()
                mesh = source_dir / (name or "inspection piece.stl")
                mesh.write_bytes(geometry(component, posed)[0])
                self.inventory.register_mesh(mesh)

            if published:
                destination = OUT / component
            elif destination is None:
                destination = scratch / "requested surveys" / f"{component} inspection"
            mesh_arg = os.path.relpath(mesh, caller) if relative else mesh
            dir_arg = os.path.relpath(destination, caller) if relative else destination
            status = run_delivered(mesh_arg, dir_arg, scratch, cwd=caller, model_only=model_only)
            report = destination / "survey.json"
            root = OUT if published else destination
            self.record(group, component, posed, regular_bytes(report, root=root), status)
            self.saved_checkpoint()
        finally:
            shutil.rmtree(scratch)

    @cache
    def exercise(self):
        assert PROGRAM.is_file(), f"Missing required program {PROGRAM}"
        initial = next(iter(self.groups["base"].values())).raw
        for group, component in LISTED_CLI.items():
            self.invoke(group, component, published=True)

        # A caller-owned destination retains each report until the next component.
        shared = scratch_dir()
        try:
            reused = shared / "current inspection"
            reused.mkdir()
            os.chown(reused, 65534, 65534)
            previous = reused / "survey.json"
            previous.write_bytes(initial or b'{"bores":[]}')
            os.chown(previous, 65534, 65534)
            self.invoke("screen-independent", "screen-leg", name=SAVED["screen-leg"],
                        model_only=True, relative=True, destination=reused)

            for group, name in zip(POSED_SCREEN, (SAVED["screen-leg"], "sample 7.stl"), strict=True):
                self.invoke(group, "screen-leg", posed=True, name=name, model_only=True)
            for component in FRESH:
                self.invoke(component, component, posed=True, model_only=True,
                            destination=reused if component == "slide" else None)
        finally:
            shutil.rmtree(shared)

        for case in REJECTED_CASES:
            self.failures[case] = rejected_input(case, initial or b'{"bores":[]}')
            self.saved_checkpoint()

        changed = self.inventory.package_fingerprint() != self.inventory.initial_package
        for component in (*SAVED, *FRESH):
            if component in REPLAYS:
                group = f"{component}-final-cli"
            elif changed:
                group = f"{component}-cli" if component in SAVED else component
            else:
                continue
            self.invoke(group, component, posed=component in FRESH,
                        listed_path=component in SAVED, model_only=component in FRESH)

    def __call__(self, key):
        if key not in SAVED:
            self.exercise()
        assert self.groups[key], f"No survey was recorded for {key}"
        return self.groups[key].values()

    @property
    def rejected(self):
        self.exercise()
        return self.failures


@pytest.fixture(scope="module")
def surveys(report_inventory):
    return Inspections(report_inventory)


@ALL_CASES
def test_survey_document(surveys, key):
    """Each requested directory must contain survey.json with the single array key bores."""
    for survey in surveys(key):
        assert set(survey.document) == {"bores"}, f"{key}: expected only the bores key"
        assert isinstance(survey.bores, list), f"{key}: expected an array of bores"


@LAND_CASES
def test_land_record_fields(surveys, key):
    """Each land object must have exactly c, d, e, r and rows, with numeric XYZ triples."""
    for survey in surveys(key):
        records = survey.bores
        assert records, f"{key}: expected records for the qualifying seats"
        for record in records:
            assert isinstance(record, dict) and set(record) == {"c", "d", "e", "r", "rows"}, (
                "A land must have exactly the five specified fields"
            )
            for name in ("c", "d", "e"):
                vector(record, name)
            radius(record)


@ALL_CASES
def test_inventory_of_qualifying_lands(surveys, key):
    """Report one land per qualifying inward cylinder, pooling its disconnected trims."""
    for survey in surveys(key):
        assert len(survey.bores) == len(survey.lands), (
            f"{key}: found {len(survey.bores)} records for {len(survey.lands)} qualifying cylinders"
        )


@LAND_CASES
def test_closest_axis_points(surveys, key):
    """The point of each matched axis closest to the input origin must be within 0.03 mm."""
    for survey in surveys(key):
        for _, land, error, _, _ in survey.pairs:
            assert error <= 0.03, f"{key}: axis near {land.c} has datum error {error:.6g} mm"


@LAND_CASES
def test_nominal_land_radii(surveys, key):
    """Each radius must be within 0.03 mm of the one-to-one matched nominal cylinder."""
    for survey in surveys(key):
        for _, land, _, error, _ in survey.pairs:
            assert error <= 0.03, f"{key}: cylinder at {land.c} has radius error {error:.6g} mm"


@LAND_CASES
def test_unoriented_axes(surveys, key):
    """Allow either sense of each bore axis within 0.03 degrees."""
    for survey in surveys(key):
        for _, land, _, _, angle in survey.pairs:
            assert angle <= 0.03, f"{key}: cylinder at {land.c} has axis error {angle:.6g} degrees"


@LAND_CASES
def test_orthonormal_clock_frames(surveys, key):
    """The reported d and e must be unit vectors and perpendicular within 0.000001."""
    for survey in surveys(key):
        records = survey.bores
        assert records, f"{key}: missing cylinder frames"
        for record in records:
            d, e = vector(record, "d"), vector(record, "e")
            errors = [abs(np.linalg.norm(v) - 1) for v in (d, e)] + [abs(d @ e)]
            assert max(errors) <= 1e-6, f"{key}: frame errors in |d|, |e| and d·e are {errors}"


@LAND_CASES
def test_station_interval_tables(surveys, key):
    """Every two-degree row contains three sorted lists of merged, positive-length intervals."""
    for survey in surveys(key):
        records = survey.bores
        assert records, f"{key}: missing interval tables"
        for record in records:
            assert isinstance(record, dict), "Each bore needs its own table of rows"
            check_rows(record.get("rows"))


@LAND_CASES
def test_access_is_on_the_reported_wall(surveys, key):
    """Plus and minus must each be a subset of the occupied wall at that clock position."""
    for survey in surveys(key):
        records = survey.bores
        assert records, f"{key}: missing access tables"
        for record in records:
            assert isinstance(record, dict), "Each bore needs an access table"
            check_access_subsets(record.get("rows"))


def assert_column(survey, column):
    errors = survey.scores[:, column]
    i = int(np.argmax(errors))
    worst = errors[i]
    land = survey.pairs[i][1]
    assert worst <= 0.015, (
        f"{survey.key}: relative interval difference {worst:.6g} at {land.c}, radius {land.r}"
    )


@LAND_CASES
def test_trimmed_wall_stations(surveys, key):
    """Wall intervals must follow the actual cylinder trims to the published 0.015 score."""
    for survey in surveys(key):
        assert_column(survey, 0)


@LAND_CASES
def test_positive_end_access(surveys, key):
    """The +d column must clear the whole shank and its leading tip on insertion from +d."""
    for survey in surveys(key):
        assert_column(survey, 1)


@LAND_CASES
def test_negative_end_access(surveys, key):
    """The −d column must clear the whole shank and its leading tip on insertion from −d."""
    for survey in surveys(key):
        assert_column(survey, 2)


REJECTED_CASES = ("permission-denied", "missing", "truncated")


@contextmanager
def watch_publication(path, original, prior_bytes):
    """Observe the open report, its logical name, and replacement of its directory."""
    libc = ctypes.CDLL(None, use_errno=True)
    fd = libc.inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
    if fd < 0:
        raise OSError(ctypes.get_errno(), "Cannot monitor the requested report")
    moved_from, moved_to, create, delete = 0x40, 0x80, 0x100, 0x200
    delete_self, move_self = 0x400, 0x800
    modify, attrib, close_write = 0x02, 0x04, 0x08
    entry_events = moved_from | moved_to | create | delete
    self_events = delete_self | move_self
    file_events = modify | attrib | close_write | self_events
    event_header = struct.Struct("iIII")
    changes = set()
    try:
        def add_watch(item, mask):
            watch = libc.inotify_add_watch(fd, os.fsencode(item), mask)
            if watch < 0:
                raise OSError(ctypes.get_errno(), f"Cannot monitor {item}")
            return watch

        file_watch = add_watch(path, file_events)
        dir_watch = add_watch(path.parent, entry_events | self_events)
        parent_watch = add_watch(path.parent.parent, entry_events)
        done, notify = os.pipe()

        def observe():
            while True:
                ready, _, _ = select([fd, done], [], [])
                try:
                    data = os.read(fd, 65536)
                except BlockingIOError:
                    if done in ready:
                        return
                    continue
                offset = 0
                while offset < len(data):
                    watch, mask, _, size = event_header.unpack_from(data, offset)
                    offset += event_header.size
                    name = data[offset:offset + size].rstrip(b"\0")
                    offset += size
                    if watch == file_watch and mask & file_events:
                        # Closing a read/write descriptor after taking a lock,
                        # or changing only metadata, need not change a report.
                        if (mask & (modify | self_events) or
                                os.pread(original.fileno(), len(prior_bytes) + 1, 0) != prior_bytes):
                            changes.add(f"report file event {mask:#x}")
                    if watch == dir_watch and mask & self_events:
                        changes.add(f"report directory event {mask:#x}")
                    if watch == dir_watch and mask & entry_events and name == os.fsencode(path.name):
                        changes.add(f"report entry event {mask:#x}")
                    if (watch == parent_watch and mask & entry_events
                            and name == os.fsencode(path.parent.name)):
                        changes.add(f"destination entry event {mask:#x}")

        try:
            with ThreadPoolExecutor(max_workers=1) as observer:
                future = observer.submit(observe)
                try:
                    yield changes
                finally:
                    os.write(notify, b"x")
                    future.result()
        finally:
            os.close(notify)
            os.close(done)
    finally:
        os.close(fd)


def rejected_input(case, prior_bytes):
    """Invoke the CLI on an unreadable model while a populated survey is published."""
    # The caller owns its working directory and may rename its destination.
    # The parent of the working directory stays fixed throughout observation.
    with tempfile.TemporaryDirectory() as anchor:
        os.chmod(anchor, 0o711)
        scratch = Path(anchor) / "caller"
        scratch.mkdir()
        os.chown(scratch, 65534, 65534)
        source = scratch / f"{case}.stl"
        if case == "permission-denied":
            shutil.copyfile(INPUTS / SAVED["screen-leg"], source)
            source.chmod(0)
        elif case == "truncated":
            release = (INPUTS / SAVED["screen-leg"]).read_bytes()
            declared = struct.unpack_from("<I", release, 80)[0]
            assert declared > 0 and len(release) == 84 + declared * FACET.itemsize, (
                "Expected a complete binary STL in the supplied print release"
            )
            source.write_bytes(release[:-FACET.itemsize])
        destination = scratch / "inspection"
        destination.mkdir()
        os.chown(destination, 65534, 65534)
        path = destination / "survey.json"
        path.write_bytes(prior_bytes)
        os.chown(path, 65534, 65534)
        with path.open("rb") as original:
            st = os.fstat(original.fileno())
            with watch_publication(path, original, prior_bytes) as changes:
                status = run_delivered(source, destination, scratch)
            return status, ((st.st_dev, st.st_ino), prior_bytes), regular_snapshot(path)[:2], changes


@pytest.mark.parametrize("case", REJECTED_CASES)
def test_unreadable_mesh_returns_failure(surveys, case):
    """An unreadable MODEL must cause a nonzero exit."""
    status, _, _, _ = surveys.rejected[case]
    assert status != 0, "The survey command accepted a model it could not read"


@pytest.mark.parametrize("case", REJECTED_CASES)
def test_unreadable_mesh_preserves_report(surveys, case):
    """An unreadable MODEL must leave the original report at its path with its contents intact."""
    _, (original_id, before), (current_id, after), changes = surveys.rejected[case]
    assert not changes, f"The published report changed during the failed invocation: {sorted(changes)}"
    assert current_id == original_id, "The command replaced or removed the previous report"
    assert after == before, "The command changed the contents of the previous report"


class InventorySnapshot(NamedTuple):
    missing: set
    extra_surveys: set
    possible_meshes: dict


class DeliveryEntry(NamedTuple):
    path: Path
    mode: int
    identity: tuple[int, int]
    content: bytes


def facet_signature(content):
    """Identify the multiset of binary STL triangles, independent of their serialization order."""
    if len(content) < 84:
        return None
    count = struct.unpack_from("<I", content, 80)[0]
    if count == 0 or len(content) != 84 + count * FACET.itemsize:
        return None
    vertices = np.frombuffer(content, FACET, count=count, offset=84)["vertices"].copy()
    if not np.isfinite(vertices).all():
        return None
    vertices[vertices == 0] = 0  # positive and negative floating-point zero are equal
    # Group the exact XYZ triples into triangles before sorting the triangles.
    points = vertices.view("V12").reshape(count, 3)
    facets = np.sort(points, axis=1).view("V36").reshape(count)
    return hashlib.sha256(np.sort(facets).tobytes()).digest()


def populated_survey(content):
    """Recognize a complete measurement report with wall data, permitting templates and maps."""
    try:
        document = json.loads(content, object_pairs_hook=unique_object)
        if not isinstance(document, dict) or set(document) != {"bores"}:
            return False
        records = document["bores"]
        if not isinstance(records, list) or not records:
            return False
        for record in records:
            if not isinstance(record, dict) or set(record) != {"c", "d", "e", "r", "rows"}:
                return False
            for key in ("c", "d", "e"):
                vector(record, key)
            radius(record)
            check_rows(record["rows"])
            if not any(row[0] for row in record["rows"]):
                return False
    except (AssertionError, ValueError, OverflowError, RecursionError):
        return False
    return True


def delivery_entries():
    """Read regular inodes without following a link or opening a device, FIFO or socket."""
    for parent, directories, names, directory_fd in os.fwalk(OUT, follow_symlinks=False):
        for name in directories + names:
            path = Path(parent) / name
            fd = os.open(name, os.O_PATH | os.O_NOFOLLOW, dir_fd=directory_fd)
            try:
                st = os.fstat(fd)
                if stat.S_ISLNK(st.st_mode):
                    content = os.fsencode(os.readlink("", dir_fd=fd))
                elif stat.S_ISREG(st.st_mode):
                    content = Path(f"/proc/self/fd/{fd}").read_bytes()
                elif stat.S_ISDIR(st.st_mode):
                    content = b""
                else:
                    continue
                yield DeliveryEntry(path, st.st_mode, (st.st_dev, st.st_ino), content)
            finally:
                os.close(fd)


class DeliveryInventory:
    def __init__(self):
        sources = [INPUTS / name for name in SAVED.values()]
        sources.extend(GRADER / "meshes" / f"{key}.stl" for key in FRESH)
        self.meshes = {facet_signature(source.read_bytes()) for source in sources}
        assert None not in self.meshes, "Cannot identify the original binary STL meshes"
        self.initial = self.snapshot()
        self.initial_package = self.package_fingerprint()

    @staticmethod
    def package_fingerprint():
        """Digest package paths, modes, link destinations and contents in stable path order."""
        entries = []
        for item in delivery_entries():
            if "__pycache__" in item.path.parts or item.path.suffix in {".pyc", ".nbc", ".nbi"}:
                continue
            relative = os.fsencode(item.path.relative_to(OUT))
            entries.append((relative, item.mode, hashlib.sha256(item.content).digest()))
        digest = hashlib.sha256()
        for relative, mode, value in sorted(entries):
            digest.update(relative + b"\0" + mode.to_bytes(4, "little") + value)
        return digest.digest()

    def register_mesh(self, mesh):
        """Include the actual posed triangles before delivering the MODEL to the command."""
        content = regular_bytes(mesh)
        assert content is not None, "Cannot fingerprint the requested MODEL"
        signature = facet_signature(content)
        assert signature is not None, "The requested MODEL is not a binary STL"
        self.meshes.add(signature)

    def snapshot(self):
        required = {PROGRAM, *(OUT / key / "survey.json" for key in SAVED)}
        missing, targets, identities = set(), set(), set()
        for path in required:
            found = regular_snapshot(path)
            if found.target is None:
                missing.add(path)
            else:
                targets.add(found.target)
                identities.add(found.identity)

        extra_surveys, possible_meshes = set(), {}
        for path, mode, identity, content in delivery_entries():
            if (not stat.S_ISREG(mode) or path in required or path in targets
                    or identity in identities):
                continue
            signature = facet_signature(content)
            if signature is not None:
                possible_meshes[path] = signature
            elif populated_survey(content):
                extra_surveys.add(path)
        return InventorySnapshot(missing, extra_surveys, possible_meshes)


@pytest.fixture(scope="module")
def report_inventory():
    return DeliveryInventory()


def test_delivery_contains_only_final_work(surveys, report_inventory):
    """Only final deliverables belong under /app/out, including resources used by the program."""
    surveys(next(iter(FINAL)))
    for phase, snapshot in (("initial", report_inventory.initial),
                            ("final", report_inventory.snapshot())):
        assert not snapshot.missing, f"Missing {phase} deliverables: {sorted(snapshot.missing)}"
        assert not snapshot.extra_surveys, f"Extra {phase} surveys: {sorted(snapshot.extra_surveys)}"
        copies = [path for path, data in snapshot.possible_meshes.items()
                  if data in report_inventory.meshes]
        assert not copies, f"Duplicated {phase} input meshes: {sorted(copies)}"
