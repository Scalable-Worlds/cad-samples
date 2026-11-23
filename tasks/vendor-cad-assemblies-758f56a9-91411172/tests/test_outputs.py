"""Acceptance measurements for trapped mould steel in the stock coordinate frame."""

import filecmp
import io
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

GRADER = Path(__file__).parent / "grader"
sys.path.insert(0, str(GRADER))
from sections import difference_lengths, solid_sections
from stock import load_stock, stock_bytes, study
from plydata import load_indexed_ply
from surfaces import buried_face, coincident_patches, surface_difference

OUT = Path("/app/out")
PYTHON = "/opt/conda/bin/python"
SETUPS = {"axial": (0, 0, 1), "side": (1, 0, 0), "inclined": (2, -1, 4)}
SOURCES = {
    "panel": Path("/app/inputs/hexapod_spherical_robot/printed_parts/upper_spherical_shell.STL"),
    "lid": Path("/app/inputs/modular_shell/node_shell_upper.STL"),
    "base": Path("/app/inputs/modular_shell/node_shell_lower.STL"),
    "lower": GRADER / "parts/lower.stl",
    "lid_b": GRADER / "parts/cover_b.stl",
    "lid_c": GRADER / "parts/cover_c.stl",
}
SOURCES.update(lid_placed=SOURCES["lid"], base_placed=SOURCES["base"])
# A second invocation uses the same stock in an independently allocated directory.
SOURCE_KEY = {"saved": "panel", **{key: key for key in SOURCES}, "repeat": "lid_c"}
T = 0.05


def scratch_dir():
    path = tempfile.mkdtemp(dir="/tmp")
    os.chown(path, 65534, 65534)
    return path


def internal_member(root, name):
    """Resolve relative links, rejecting a cycle or any hop outside the package."""
    def walk(here, parts, active):
        for item in parts:
            if not here.is_dir():
                return None
            if item == "..":
                if here == root:
                    return None
                here = here.parent
                continue
            here = here / item
            if here.is_symlink():
                link = Path(os.readlink(here))
                if link.is_absolute() or here in active:
                    return None
                here = walk(here.parent, link.parts, active | {here})
                if here is None:
                    return None
        return here

    member = Path(name)
    if root.is_symlink() or member.is_absolute():
        return None
    return walk(root, member.parts, set())


def regular_member(root, name):
    path = internal_member(root, name)
    return path if path is not None and path.is_file() else None


def run_delivered(command, scratch, destination):
    """Execute one invocation, then retain only safe locations of stopped writers' output."""
    members = {}
    try:
        status = subprocess.run(
            [sys.executable, "-I", str(GRADER / "invoke.py"), scratch, *command],
            stdin=subprocess.DEVNULL,
        ).returncode
    finally:
        # The dedicated launcher has already reaped the entire descendant tree.
        # Record in-directory targets, then remove links before opening any files.
        for setup in SETUPS:
            name = setup + ".ply"
            members[name] = regular_member(destination, name)
        subprocess.run(["find", scratch, str(OUT), "-type", "l", "-delete"], check=True)
    return status, members


def output_tail(scratch, limit=400):
    """Read the last characters after subprocess cleanup."""
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


# Additional CAD exports are study products. Implementation modules and their
# resources can accompany the command in the saved submission.
CAD_SUFFIXES = {".ply", ".stl", ".obj", ".off", ".3mf", ".amf", ".step", ".stp",
                ".iges", ".igs", ".brep", ".fcstd", ".gltf", ".glb", ".dxf"}


def capture_study(directory, archive, members=None):
    archive.mkdir()
    if members is None:
        members = {f"{s}.ply": regular_member(directory, f"{s}.ply") for s in SETUPS}
    expected = set(members.values()) | {directory / k for k in members}
    extras = []
    if directory.is_dir() and not directory.is_symlink():
        extras = sorted(str(p.relative_to(directory)) for p in directory.rglob("*")
                        if p.suffix.lower() in CAD_SUFFIXES and p not in expected and p.is_file())
    for name, path in members.items():
        if path is not None and path.is_file():
            shutil.copyfile(path, archive / name)
    return extras


@contextmanager
def invocation_package():
    """Keep the submission private while the published path holds a disposable copy."""
    with tempfile.TemporaryDirectory(dir=OUT.parent) as private:
        package = Path(private) / "package"
        package.mkdir()
        entries = list(OUT.iterdir())
        moved = []
        live = False
        mode = OUT.stat().st_mode
        try:
            for entry in entries:
                shutil.move(entry, package / entry.name)
                moved.append(entry.name)
            live = True
            shutil.copytree(package, OUT, dirs_exist_ok=True, symlinks=True)
            os.chmod(OUT, mode)
            for path in (OUT, *OUT.rglob("*")):
                if path.is_symlink() and internal_member(OUT, path.relative_to(OUT)) is None:
                    path.unlink()
                else:
                    os.chown(path, 65534, 65534, follow_symlinks=False)
            yield
        finally:
            if live:
                for child in OUT.iterdir():
                    if child.is_symlink() or not child.is_dir():
                        child.unlink()
                    else:
                        shutil.rmtree(child)
            for name in moved:
                shutil.move(package / name, OUT / name)
            os.chmod(OUT, mode)


@pytest.fixture(scope="module")
def results():
    with tempfile.TemporaryDirectory() as private:
        private = Path(private)
        saved_extras = capture_study(OUT, private / "saved")

        @lru_cache(None)
        def obtain(case):
            archive = private / case
            if case == "saved":
                return archive, 0, "", saved_extras
            program = regular_member(OUT, "pull.py")
            if program is None:
                return archive, None, "missing /app/out/pull.py", []
            scratch = scratch_dir()
            try:
                working = Path(scratch)
                part, destination = working / "stock.STL", working / "steel"
                key = SOURCE_KEY[case]
                part.write_bytes(stock_bytes(SOURCES[key], key))
                os.chmod(part, 0o644)
                destination.mkdir()
                os.chown(destination, 65534, 65534)
                with invocation_package():
                    status, members = run_delivered(
                        [PYTHON, str(OUT / "pull.py"), "--part", str(part), "--out", str(destination)],
                        scratch, destination,
                    )
                    if status:
                        return archive, status, output_tail(scratch), []
                    extras = capture_study(destination, archive, members)
                    return archive, 0, "", extras
            finally:
                shutil.rmtree(scratch)

        yield obtain


def require_result(results, case):
    directory, status, detail, _ = results(case)
    assert status == 0, f"Cannot generate {case} study (exit {status}). {detail}"
    return directory


@pytest.fixture(scope="module", params=tuple(SOURCE_KEY))
def case(request):
    return request.param


@pytest.fixture(scope="module", params=tuple(SETUPS))
def setup(request):
    return request.param


@lru_cache(maxsize=2)
def read_ply(path):
    assert path.is_file(), f"Missing triangle PLY at {path.name}"
    try:
        return load_indexed_ply(path)
    except (ValueError, TypeError, KeyError, OverflowError, OSError) as error:
        assert False, f"Cannot parse {path.name} as indexed triangle PLY: {error}"


@pytest.fixture(scope="module")
def mesh(results, case, setup):
    return read_ply(require_result(results, case) / (setup + ".ply"))


@lru_cache(None)
def stock(key):
    return load_stock(io.BytesIO(stock_bytes(SOURCES[key], key)))


@lru_cache(None)
def nominal(key, setup):
    return study(stock(key), SETUPS[setup])


def project(mesh, basis):
    return (np.asarray(mesh.vertices) @ basis)[mesh.faces]


@lru_cache(maxsize=3)
def nominal_lines(key, setup):
    spec = nominal(key, setup)
    points = np.vstack((spec.lattice, spec.cell_points))
    return solid_sections(project(stock(key), spec.basis), points, bounded_air=True)


@lru_cache(maxsize=3)
def delivered_lines(path, key, setup):
    mesh = read_ply(path)
    spec = nominal(key, setup)
    return solid_sections(project(mesh, spec.basis), np.vstack((spec.lattice, spec.cell_points)))


@pytest.fixture(scope="module")
def deviations(results, case, setup):
    path = require_result(results, case) / (setup + ".ply")
    key = SOURCE_KEY[case]
    actual = delivered_lines(path, key, setup)
    expected = nominal_lines(key, setup)
    distances = difference_lengths(expected, actual)
    split = len(nominal(key, setup).lattice)
    return distances[:split], distances[split:]


def metric(name, value, bar, sense="max"):
    print(f"METRIC name={name} value={float(value):.16g} bar={float(bar):.16g} sense={sense}")


def test_final_study_inventory(results, case):
    """“Write only your final deliverables” and “write all three files there.”"""
    directory = require_result(results, case)
    assert all((directory / f"{s}.ply").is_file() for s in SETUPS), "Missing a required study"
    extras = results(case)[3]
    assert not extras, f"Additional CAD products in the output destination: {extras}"


def test_union_boundary_has_no_internal_interfaces(mesh):
    """“Deliver its boundary” excludes positive-area interfaces between joined steel pieces."""
    overlap = coincident_patches(mesh)
    assert overlap is None, f"Faces {overlap} overlap in area or cross through each other’s interiors"


def test_indexed_triangle_ply(mesh):
    """“Deliver its boundary as an indexed triangle PLY” in each named setup."""
    assert mesh.faces.ndim == 2 and mesh.faces.shape[1] == 3, "Expected a table of indexed triangles"


def test_empty_or_required_faces(mesh, case, setup):
    """“Zero trapped volume requires a PLY with zero faces. Above 0.1 mm³ the PLY must contain faces.”"""
    volume = nominal(SOURCE_KEY[case], setup).volume
    if volume == 0:
        assert len(mesh.faces) == 0, "Stock has no trapped steel for this pull"
    elif volume > 0.1:
        assert len(mesh.faces) > 0, "Required trapped-steel boundary has no faces"


def test_indexed_watertightness(mesh):
    """“For each nonempty PLY, require watertightness” using stored vertex indices."""
    assert not len(mesh.faces) or mesh.is_watertight, "Open or non-manifold indexed boundary"


def test_adjacent_face_winding(mesh):
    """“Consistent outward winding” requires opposite traversal across shared edges."""
    assert not len(mesh.faces) or mesh.is_winding_consistent, "Inconsistent orientation across an indexed edge"


def test_positive_signed_volume(mesh, case, setup):
    """A nonempty PLY needs “positive signed volume using the stored vertex indices.”"""
    if len(mesh.faces):
        metric(f"{case}_{setup}_signed_mm3", mesh.volume, 0, "min")
        assert mesh.volume > 0, f"Nonpositive signed volume {mesh.volume}"


def test_outward_winding_of_each_boundary(mesh):
    """“Consistent outward winding” places air on the outward side, including inside cavities."""
    error = buried_face(mesh)
    assert error is None, error


def test_stock_bounds(mesh, case, setup):
    """“All vertices must be inside the stock's XYZ bounds expanded by 0.05 mm.”"""
    bounds = stock(SOURCE_KEY[case]).bounds
    vertices = np.asarray(mesh.vertices)
    excess = max(np.max(bounds[0] - vertices, initial=0), np.max(vertices - bounds[1], initial=0))
    metric(f"{case}_{setup}_bounds_mm", excess, T)
    assert excess <= T, f"Vertex outside expanded stock bounds by {excess - T:.8g} mm"


def test_integrated_gap_volume(mesh, case, setup):
    """“Absolute volume error must be at most 0.05 mm × A.”"""
    expected = nominal(SOURCE_KEY[case], setup)
    error = abs((mesh.volume if len(mesh.faces) else 0.0) - expected.volume)
    bar = T * expected.rectangle_area
    metric(f"{case}_{setup}_volume_mm3", error, bar)
    assert error <= bar, f"Trapped volume error {error:.8g} mm³ exceeds {bar:.8g} mm³"


def test_lattice_sections(deviations, case, setup):
    """Use the “square lattice of 0.5 mm pitch” for the 0.05 mm line-deviation requirement."""
    maximum = np.max(deviations[0], initial=0.0)
    metric(f"{case}_{setup}_lattice_mm", maximum, T)
    assert maximum <= T, f"Lattice line {np.argmax(deviations[0])} differs by {maximum:.8g} mm"


def test_projected_cell_sections(deviations, case, setup):
    """“Also measure at shapely.point_on_surface of each polygonized cell” with the stated clearance."""
    maximum = np.max(deviations[1], initial=0.0)
    metric(f"{case}_{setup}_cells_mm", maximum, T)
    assert maximum <= T, f"Projected cell {np.argmax(deviations[1])} differs by {maximum:.8g} mm"


@pytest.mark.parametrize("original,again", [("saved", "panel"), ("lid_c", "repeat")],
                         ids=["saved-study", "fresh-destinations"])
def test_reproducible_study(results, setup, original, again):
    """“Generate it with” the stated command, which must write all three files “deterministically.”"""
    a, b = [require_result(results, case) / (setup + ".ply") for case in (original, again)]
    assert a.is_file() and b.is_file(), f"Cannot compare missing {setup} studies"
    if filecmp.cmp(a, b, shallow=False):
        return
    error = surface_difference(read_ply(a), read_ply(b))
    assert error is None, error
