"""Compare laser and brake geometry in the datum frame of a registered folded master."""

import filecmp
import glob
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from functools import cached_property, lru_cache

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import ezdxf
import numpy as np
import pytest
import shapely
from scipy.optimize import linear_sum_assignment

GRADER = Path(__file__).parent / "grader"
OUT = Path("/app/out")
INPUTS = Path("/app/inputs")
PYTHON = "/opt/conda/bin/python"
TOLERANCE = 0.15
UNPRIVILEGED = dict(user=65534, group=65534, extra_groups=[])
JOBS = GRADER / "jobs"
CASES = ["saved_base", "regenerated_base",
         *(p.stem for p in sorted(JOBS.glob("*.json")) if p.stem != "saved_base")]


def scratch_dir():
    path = tempfile.mkdtemp(dir="/tmp")
    os.chown(path, 65534, 65534)
    return path


def still_running():
    """Whether a process of UID 65534 is alive (a zombie is not)."""
    for path in glob.glob("/proc/[0-9]*/task/[0-9]*/status"):
        try:
            with open(path) as file:
                status = file.read()
        except OSError:
            continue
        if "\nUid:\t65534\t" in status and "\nState:\tZ" not in status:
            return True
    return False


def release_target(scratch, destination):
    """Resolve the output as its author after exit, without opening its contents."""
    result = subprocess.run(
        ["/usr/bin/realpath", "-ze", "--", destination, destination / "blank.dxf"],
        cwd="/", stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        **UNPRIVILEGED,
    )
    if result.returncode == 0:
        release, drawing, _ = result.stdout.split(b"\0")
        release = Path(os.fsdecode(release))
        drawing = Path(os.fsdecode(drawing))
        if release.is_relative_to(scratch) and drawing.is_relative_to(release):
            return drawing
    return None


@contextmanager
def private_calibration_stock():
    """Keep the calibration files private while the CLI has its designated STL copy."""
    mode = INPUTS.stat().st_mode & 0o7777
    INPUTS.chmod(0o700)
    try:
        yield
    finally:
        INPUTS.chmod(mode)


def run_delivered(command, scratch, destination):
    """Run the CLI with only its scratch stock and stop writers before resolving the DXF."""
    env = {**os.environ, "HOME": scratch, "TMPDIR": scratch, "MPLCONFIGDIR": scratch}
    env.pop("PYTEST_CURRENT_TEST", None)
    with private_calibration_stock(), open(f"{scratch}/stdout", "wb") as out, open(f"{scratch}/stderr", "wb") as err:
        try:
            result = subprocess.run(
                command, cwd=scratch, env=env, stdin=subprocess.DEVNULL, stdout=out,
                stderr=err, **UNPRIVILEGED,
            ).returncode
        finally:
            subprocess.run(["sh", "-c", "kill -9 -1"], env=env, **UNPRIVILEGED)
            while still_running():
                time.sleep(0.01)
            try:
                target = release_target(scratch, destination)
            finally:
                subprocess.run(["find", scratch, "/app/out", "-type", "l", "-delete"], check=True)
    return result, target


def output_tail(scratch, limit=400):
    """Read the end of the delivered program's log after cleanup."""
    if limit <= 0:
        return ""
    for name in ("stderr", "stdout"):
        path = f"{scratch}/{name}"
        if not os.path.isfile(path):
            continue
        with open(path, "rb") as log:
            end = log.seek(0, os.SEEK_END)
            if end:
                log.seek(max(0, end - 4 * limit))
                return log.read(4 * limit).decode("utf-8", errors="replace")[-limit:]
    return ""


def read_drawing(path):
    assert path.is_file(), f"Missing drawing {path.name}"
    try:
        return ezdxf.readfile(path)
    except (OSError, ValueError, ezdxf.DXFError) as exc:
        assert False, f"Cannot reopen {path.name} as DXF: {exc}"


def generate_drawing(stl, datum, existing=False):
    script = OUT / "cutblank.py"
    assert script.is_file(), "Missing cutblank.py"
    scratch = scratch_dir()
    try:
        local_part = Path(scratch) / "formed stock.STL"
        shutil.copyfile(stl, local_part)
        destination = Path(scratch) / "release"
        if existing:
            destination.mkdir()
            os.chown(destination, 65534, 65534)
        command = [PYTHON, str(script), "--part", str(local_part), "--datum",
                   *map(str, datum), "--out", str(destination)]
        result, target = run_delivered(command, scratch, destination)
        assert result == 0, f"Blank generation exited {result}. {output_tail(scratch)}"
        assert target is not None, "Missing blank.dxf within the requested release directory"
        return read_drawing(target)
    finally:
        shutil.rmtree(scratch)


@dataclass
class PlacedEntity:
    source: object
    matrix: object
    roundoff: object = None

    @cached_property
    def local_vertices(self):
        if self.source.dxftype() == "LINE":
            return (self.source.dxf.start, self.source.dxf.end)
        if self.source.dxftype() == "ATTRIB":
            return (self.source.ocs().to_wcs(self.source.dxf.insert),)
        return tuple(self.source.vertices_in_wcs())

    @cached_property
    def vertices(self):
        """Get path vertices or the label insertion point in modelspace."""
        if self.matrix is None:
            return self.local_vertices
        return tuple(self.matrix.transform_vertices(self.local_vertices))

    def in_xy_plane(self):
        if all(p.z == 0 for p in self.vertices):
            return True
        if self.matrix is None:
            return False
        # Account for the operations that actually placed these vertices. A
        # half-turn, for example, need not evaluate to a bitwise-zero sine.
        return all(abs(world.z) <= np.dot([abs(v) for v in local] + [1], self.roundoff)
                   for world, local in zip(self.vertices, self.local_vertices))


def visible_layer(doc, name):
    if name not in doc.layers:
        return True
    layer = doc.layers.get(name)
    return not (layer.is_off() or layer.is_frozen())


def effective_entities(entities, inherited="0", ancestry=(), parent_matrix=None,
                       magnitude=None, visible=True):
    """Expand modelspace block references, including layer-0 inheritance."""
    for entity in entities:
        layer = entity.dxf.layer.upper()
        if layer == "0":
            layer = inherited
        if entity.dxftype() != "INSERT":
            if not ancestry or entity.dxftype() != "ATTDEF":
                roundoff = (32 * np.finfo(float).eps * len(ancestry) * magnitude
                            if ancestry else None)
                yield layer, PlacedEntity(entity, parent_matrix, roundoff)
            continue
        name = entity.dxf.name.casefold()
        assert name not in ancestry, "Cyclic block reference in drawing"
        inserts = entity.multi_insert() if entity.mcount > 1 else [entity]
        for insert in inserts:
            block = insert.block()
            assert block is not None, f"Missing block definition for {insert.dxf.name}"
            shown = visible and not insert.dxf.invisible and visible_layer(insert.doc, layer)
            mode = insert.doc.header.get("$ATTMODE", 1)
            for attr in insert.attribs:
                attr_layer = attr.dxf.layer.upper()
                if attr_layer == "0":
                    attr_layer = layer
                if (shown and mode != 0 and (mode == 2 or not attr.is_invisible)
                        and not attr.dxf.invisible and visible_layer(insert.doc, attr_layer)):
                    # Attached labels already have coordinates in the INSERT's
                    # owning layout. Only its ancestors still act on them.
                    yield attr_layer, PlacedEntity(attr, parent_matrix)
            # A rotation inside a nonuniform scale can produce a shear. Keep the
            # full transform until the final coordinates are needed.
            matrix = insert.matrix44()
            # Sum absolute affine terms before composition so cancellations,
            # nested scales and translations are represented in the error bound.
            weights = np.abs(np.asarray(list(matrix.rows()))) @ (
                np.ones(4) if magnitude is None else magnitude)
            if parent_matrix is not None:
                matrix = matrix @ parent_matrix
            yield from effective_entities(block, layer, (*ancestry, name), matrix, weights, shown)


@lru_cache
def master_developer():
    spec = importlib.util.spec_from_file_location("nominal", GRADER / "nominal.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.develop


def closed_coords(placed):
    vertices = np.array([tuple(p) for p in placed.vertices], dtype=float)
    assert len(vertices) >= 3, "CUT contains a loop with fewer than three vertices"
    assert np.all(np.isfinite(vertices)), "Non-finite CUT coordinate"
    assert placed.source.closed or np.array_equal(vertices[0], vertices[-1]), "Open CUT polyline"
    return vertices[:, :2]


class Contours:
    def __init__(self, coords):
        self.rings = [shapely.LinearRing(xy) for xy in coords]
        self.polygons = [shapely.Polygon(r) for r in self.rings]

    @cached_property
    def depths(self):
        # Whole-polygon containment is insensitive to winding and starting vertex.
        return [sum(other.contains(p) for j, other in enumerate(self.polygons) if i != j)
                for i, p in enumerate(self.polygons)]

    def at_depth(self, depth):
        return [r for r, d in zip(self.rings, self.depths) if d == depth]


@dataclass
class Drawing:
    key: str
    doc: object
    expected: dict

    @cached_property
    def layers(self):
        grouped = {layer: [] for layer in ("CUT", "INNER", "OUTER")}
        for layer, entity in effective_entities(self.doc.modelspace()):
            if layer in grouped:
                grouped[layer].append(entity)
        return grouped

    @property
    def cuts(self):
        entities = self.layers["CUT"]
        assert entities, "Missing CUT geometry in modelspace"
        assert all(e.source.dxftype() == "LWPOLYLINE" for e in entities), "CUT requires LWPOLYLINEs"
        return entities

    @property
    def brakes(self):
        entities = self.layers["INNER"] + self.layers["OUTER"]
        assert all(e.source.dxftype() == "LINE" for e in entities), "Brake layers require LINEs"
        return entities

    @cached_property
    def contours(self):
        return Contours([closed_coords(e) for e in self.cuts])

    @cached_property
    def nominal_contours(self):
        return Contours(self.expected["CUT"])


@pytest.fixture(scope="module", params=CASES)
def drawing(request):
    name = request.param
    geometry = "saved_base" if name == "regenerated_base" else name
    job = json.loads((JOBS / f"{geometry}.json").read_text())
    stl = INPUTS / job["part"]
    if not stl.is_file():
        stl = GRADER / "stock" / job["part"]
    expected = master_developer()(stl.with_suffix(".STEP"), stl, job["datum"])
    if name == "saved_base":
        doc = read_drawing(OUT / "blank.dxf")
    else:
        doc = generate_drawing(stl, job["datum"], existing=(name == "regenerated_base"))
    return Drawing(name, doc, expected)


def bottleneck(cost):
    """Find the bijection minimizing its largest entry, including tied costs."""
    levels = np.unique(cost)
    lo, hi = 0, len(levels) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        blocked = cost > levels[mid]
        rows, cols = linear_sum_assignment(blocked)
        if blocked[rows, cols].any():
            lo = mid + 1
        else:
            hi = mid
    rows, cols = linear_sum_assignment(cost > levels[lo])
    return levels[lo], rows, cols


def report_limit(name, value, what):
    print(f"METRIC name={name} value={value} bar={TOLERANCE} sense=max", flush=True)
    assert value <= TOLERANCE, f"{what}: {value:.6g} mm exceeds {TOLERANCE} mm"


def compare_curves(actual, nominal, name):
    assert len(actual) == len(nominal), f"{name}: expected {len(nominal)} curves, got {len(actual)}"
    if not nominal:
        return
    # Extrema bound Hausdorff distance from below, even with discrete densification.
    # One feasible assignment supplies an upper bound. Only pairs that can improve
    # on it need the expensive symmetric Hausdorff calculation.
    a = shapely.bounds(actual)
    b = shapely.bounds(nominal)
    bound = np.abs(a[:, None, :] - b[None, :, :]).max(axis=2)
    assert np.isfinite(bound).all(), "Non-finite curve bounds"
    lower, rows, cols = bottleneck(bound)
    if lower > TOLERANCE:
        report_limit(name + "_extent_lower_bound", lower, "Best possible pairing, extent lower bound")
    cost = np.full_like(bound, np.inf)

    def measure(i, j):
        return shapely.hausdorff_distance(actual[i], nominal[j], densify=0.25)

    for i, j in zip(rows, cols):
        cost[i, j] = measure(i, j)
    upper = cost[rows, cols].max()
    for i, j in zip(*np.where((bound <= upper) & np.isinf(cost))):
        cost[i, j] = measure(i, j)
    error, _, _ = bottleneck(cost)
    report_limit(name, error, "Maximum Hausdorff distance in the best one-to-one pairing")


def test_final_release_inventory():
    """'Only your final deliverables' excludes duplicate folded inputs and surplus copies of the release."""
    saved = OUT / "blank.dxf"
    assert saved.is_file() and (OUT / "cutblank.py").is_file(), "Missing final script or drawing"
    sources = list(INPUTS.iterdir())
    for path in OUT.rglob("*"):
        if path == saved or not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix in (".stl", ".step", ".stp"):
            assert not any(filecmp.cmp(path, src, shallow=False) for src in sources), (
                f"Surplus copy of folded stock: {path.relative_to(OUT)}"
            )
        if suffix == ".dxf":
            assert not filecmp.cmp(path, saved, shallow=False), f"Duplicate release drawing: {path.relative_to(OUT)}"


def test_modelspace_units(drawing):
    """'DXF modelspace XY in mm at 1:1' requires compatible drawing units."""
    assert drawing.layers["CUT"], "No modelspace CUT release"
    assert drawing.doc.header.get("$INSUNITS", 0) in (0, 4), "DXF declares units other than mm"


def test_one_polyline_per_cut_loop(drawing):
    """'One ... LWPOLYLINE per loop' covers the perimeter and every hole."""
    assert len(drawing.cuts) == len(drawing.expected["CUT"]), "Incorrect number of CUT polylines"


def test_closed_cut_paths(drawing):
    """'One simple closed ... LWPOLYLINE per loop' requires geometric closure."""
    for entity in drawing.cuts:
        closed_coords(entity)


def test_straight_cut_segments(drawing):
    """'Straight-segment LWPOLYLINE' excludes circular bulge segments."""
    for cut in drawing.cuts:
        points = cut.source.get_points("b")
        active = points if cut.source.closed else points[:-1]
        assert all(b == 0 for b, in active), "CUT path includes a curved segment"


def test_simple_cut_loops(drawing):
    """'Simple closed' cut loops must enclose area without self-intersection."""
    for ring, polygon in zip(drawing.contours.rings, drawing.contours.polygons):
        assert ring.is_simple and polygon.is_valid and polygon.area > 0, "Self-crossing or degenerate CUT loop"


def test_disjoint_cut_boundaries_and_nesting(drawing):
    """'Without duplicate cuts' and 'preserving nesting' require distinct nested loops."""
    contours = drawing.contours
    for i, ring in enumerate(contours.rings):
        assert all(ring.disjoint(other) for other in contours.rings[:i]), "CUT boundaries cross, touch or overlap"
    assert sorted(contours.depths) == sorted(drawing.nominal_contours.depths), "Incorrect cut-loop nesting"


def test_neutral_blank_perimeter(drawing):
    """The 'perimeter' must match the nominal mid-gauge development in drawing coordinates."""
    compare_curves(drawing.contours.at_depth(0), drawing.nominal_contours.at_depth(0),
                   drawing.key + "_perimeter_mm")


def test_datum_registered_hole_profiles(drawing):
    """'Every hole' must match its nominal curve in the datum-defined axes and origin."""
    actual, nominal = drawing.contours, drawing.nominal_contours
    for depth in sorted(set(actual.depths + nominal.depths) - {0}):
        compare_curves(actual.at_depth(depth), nominal.at_depth(depth), f"{drawing.key}_holes_{depth}_mm")


def test_brake_line_entities(drawing):
    """'One LINE per uninterrupted positive-length interval' requires straight line entities."""
    for line in drawing.brakes:
        start, end = np.asarray(line.vertices)
        assert np.all(np.isfinite([start, end])) and np.any(start != end), "Invalid brake LINE"


def test_signed_material_intervals(drawing):
    """'INNER ... concave' and 'OUTER ... convex' distinguish the picked skin at each bend."""
    for layer in ("INNER", "OUTER"):
        assert len(drawing.layers[layer]) == len(drawing.expected[layer]), f"Incorrect {layer} interval count"


def test_clipped_half_angle_generators(drawing):
    """The 'generator halfway in angle ... clipped to material' fixes each signed interval."""
    for layer in ("INNER", "OUTER"):
        actual = []
        for e in drawing.layers[layer]:
            assert e.source.dxftype() == "LINE", f"{layer} contains a non-LINE mark"
            actual.append(shapely.LineString([p.xy for p in e.vertices]))
        nominal = [shapely.LineString(pair) for pair in drawing.expected[layer]]
        compare_curves(actual, nominal, f"{drawing.key}_{layer.lower()}_mm")


def test_zero_elevation(drawing):
    """'Modelspace XY ... with zero elevation' places the entire release in the drawing plane."""
    for placed in drawing.cuts + drawing.brakes:
        assert placed.in_xy_plane(), "Release geometry is above or below XY"
        assert placed.source.dxf.get("thickness", 0) == 0, "Extruded drawing entity"


def test_zero_polyline_widths(drawing):
    """'Zero widths' requires zero geometric width along every active CUT segment."""
    for cut in drawing.cuts:
        entity = cut.source
        assert entity.dxf.get("const_width", 0) == 0, "Nonzero constant CUT width"
        widths = entity.get_points("se")
        active = widths if entity.closed else widths[:-1]
        assert all(s == 0 and e == 0 for s, e in active), "Nonzero segment width"
