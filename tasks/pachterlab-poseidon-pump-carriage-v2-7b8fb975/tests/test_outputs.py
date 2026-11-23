"""Verify print-fit geometry against the original, trimmed IGES stock."""

import os
from functools import cache, cached_property
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import FreeCAD
import Part
import pytest

DATA = Path(__file__).resolve().parent / "grader"
sys.path.insert(0, str(DATA))
from supports import groups_for, load_original
from occupancy import audit

OUT = Path("/app/out")
INPUT = Path("/app/inputs")
SOURCES = {
    "carriage": INPUT / "Pump Carriage.iges",
    "brace": INPUT / "Pump Syringe Brace.iges",
    "slide": DATA / "slide.iges",
    "mount": DATA / "mount.iges",
    "scope_base": DATA / "scope_base.iges",
    "split_brace": DATA / "brace_sectors.iges",
    "stage": DATA / "stage.iges",
    "carriage_rerun": INPUT / "Pump Carriage.iges",
}
ALL_PARTS = pytest.mark.parametrize("part", tuple(SOURCES), indirect=True, scope="module")
WITH_SEATS = pytest.mark.parametrize(
    "part", [name for name in SOURCES if name != "slide"], indirect=True, scope="module"
)
CAD_SUFFIXES = {".step", ".stp", ".iges", ".igs", ".brep", ".brp", ".stl",
                ".fcstd", ".3mf", ".x_t", ".x_b", ".sat", ".sab", ".dxf", ".dwg",
                ".sldprt", ".sldasm", ".catpart", ".catproduct", ".stpz", ".stepz"}


def tree_files(root):
    """Visit entry paths without descending through directory links."""
    for folder, directories, names in os.walk(root, followlinks=False):
        base = Path(folder)
        for name in names:
            yield base / name
        for name in directories:
            if (base / name).is_symlink():
                yield base / name


def local_file(path, root):
    """Return a regular target inside the specified real directory."""
    try:
        target = path.resolve(strict=True)
        if target.is_relative_to(root) and target.is_file():
            return target
    except (OSError, RuntimeError):
        pass
    return None


def copy_submission(source, destination):
    destination.mkdir()
    for folder, _, names in os.walk(source, followlinks=False):
        dest = destination / Path(folder).relative_to(source)
        dest.mkdir(parents=True, exist_ok=True)
        for name in names:
            target = local_file(Path(folder) / name, source)
            if target is not None:
                shutil.copy2(target, dest / name)


def cad_product(name):
    path = Path(name)
    if path.suffix.lower() in {".gz", ".bz2", ".xz"}:
        path = Path(path.stem)
    return path.suffix.lower() in CAD_SUFFIXES


def extra_models(directory, primary=None):
    result = directory / "fit.step"
    primary = primary or local_file(result, directory)
    return sorted(str(path.relative_to(directory)) for path in tree_files(directory)
                  if cad_product(path) and path != result
                  and (primary is None or local_file(path, directory) != primary))


def log_tail(stream):
    """Return recent diagnostics from the original open file descriptor."""
    size = stream.seek(0, os.SEEK_END)
    stream.seek(max(0, size - 16384))
    return "\n".join(stream.read().decode("utf-8", errors="replace").splitlines()[-20:])


def editor_user_is_active():
    for status_path in Path("/proc").glob("[0-9]*/task/[0-9]*/status"):
        try:
            status = status_path.read_bytes()
        except OSError:
            continue
        if b"\nUid:\t65534\t" in status and b"\nState:\tZ" not in status:
            return True
    return False


def run_editor(submission, source, destination, layout):
    assert (submission / "relieve.py").is_file(), "Missing /app/out/relieve.py"
    with tempfile.TemporaryDirectory(dir="/tmp") as directory:
        # Keep the scratch pathname stable even when its contents are writable.
        parent = Path(directory)
        parent.chmod(0o711)
        scratch = parent / "work"
        scratch.mkdir()
        os.chown(scratch, 65534, 65534)
        implementation = scratch / ("out" if layout else "editor")
        def omit_model_files(folder, names):
            return [name for name in names if cad_product(name) and (Path(folder) / name).is_file()]

        shutil.copytree(submission, implementation, ignore=omit_model_files)
        for item in [implementation, *implementation.rglob("*")]:
            os.chown(item, 65534, 65534)
        inputs = scratch / "inputs"
        inputs.mkdir()
        source_copy = inputs / source.name
        shutil.copyfile(source, source_copy)
        out = implementation if layout else scratch / "revised"
        if not layout:
            out.mkdir()
            os.chown(out, 65534, 65534)
        command = ["/opt/conda/bin/python", str(implementation / "relieve.py"),
                   "--part", str(source_copy), "--out", str(out)]
        env = {**os.environ, "HOME": str(scratch), "TMPDIR": str(scratch),
               "MPLCONFIGDIR": str(scratch)}
        env.pop("PYTEST_CURRENT_TEST", None)
        # The saved STEP, development masters and any other release files stay
        # in the private grading tree while the chosen input is being edited.
        app = OUT.parent
        old_mode = stat.S_IMODE(app.stat().st_mode)
        with (parent / "stdout").open("w+b") as stdout, (parent / "stderr").open("w+b") as stderr:
            try:
                app.chmod(0o700)
                rc = subprocess.run([sys.executable, "-I", str(DATA / "invoke.py"), *command],
                                    cwd=scratch, env=env, stdin=subprocess.DEVNULL,
                                    stdout=stdout, stderr=stderr).returncode
            finally:
                try:
                    subprocess.run(["sh", "-c", "kill -9 -1"], env=env,
                                   user=65534, group=65534, extra_groups=[])
                    while editor_user_is_active():
                        time.sleep(0.01)
                    # Resolve a local entry point using metadata only. Read the
                    # regular file after cleanup, even when the entry was a link.
                    target = None if out.is_symlink() else local_file(out / "fit.step", out)
                    subprocess.run(["find", scratch, OUT, "-type", "l", "-delete"], check=True)
                finally:
                    app.chmod(old_mode)
            assert rc == 0, (
                f"Editor failed on {source.name} (exit {rc})\n"
                f"stdout (last lines)\n{log_tail(stdout)}\n"
                f"stderr (last lines)\n{log_tail(stderr)}"
            )
        assert out.is_dir(), "The destination directory was replaced"
        assert target is not None, "The requested destination has no local, regular fit.step"
        extras = extra_models(out, target)
        shutil.copyfile(target, destination)
        return extras


class Inspection:
    def __init__(self, source, step, extras):
        self.source = source
        self.extras = extras
        assert step.is_file(), "Missing /app/out/fit.step"
        try:
            self.export = Part.read(str(step))
        except Exception as exc:
            raise AssertionError(f"FreeCAD Part.read could not reopen fit.step: {exc}") from exc
        assert not self.export.isNull(), "STEP import returned a null shape"

    @cached_property
    def solid(self):
        solids = self.export.Solids
        assert len(solids) == 1, f"Expected one solid in the STEP, found {len(solids)}"
        return solids[0]

    @cached_property
    def occupancy(self):
        nominal = load_original(self.source)
        return audit(nominal, groups_for(nominal), self.solid)


@pytest.fixture(scope="module")
def inspect_part():
    with tempfile.TemporaryDirectory(dir="/tmp") as directory:
        private = Path(directory)
        submission = private / "submitted"
        copy_submission(OUT, submission)

        @cache
        def inspect(name):
            source = SOURCES[name]
            if name == "carriage":
                step = submission / "fit.step"
                extras = extra_models(OUT)
            else:
                case_dir = private / name
                case_dir.mkdir(exist_ok=True)
                step = case_dir / "fit.step"
                extras = run_editor(submission, source, step, name == "carriage_rerun")
            return Inspection(source, step, extras)

        yield inspect


@pytest.fixture(scope="module")
def part(request, inspect_part):
    return inspect_part(request.param)


def root_components(shape):
    if shape.ShapeType in ("Compound", "CompSolid"):
        for child in shape.childShapes():
            yield from root_components(child)
    else:
        yield shape.ShapeType


def check_occupancy(part, key):
    result = part.occupancy[key]
    assert result["tested"] > 0, f"No applicable {key} material measurements"
    assert result["errors"] == 0, (
        f'{result["errors"]} / {result["tested"]} mismatches for {key}. '
        f'First disagreement: {result["first"]}'
    )


@ALL_PARTS
def test_one_material_solid(part):
    """“Export the full resulting material as one … solid in a STEP file.”"""
    components = list(root_components(part.export))
    assert components == ["Solid"], f"Expected a single material solid, found {components}"


@ALL_PARTS
def test_valid_brep(part):
    """“Reopen the result with FreeCAD Part.read, check validity.”"""
    assert part.export.isValid(), "The exchanged B-rep is invalid"


@ALL_PARTS
def test_closed_shells(part):
    """“One valid, closed, positive-volume B-rep solid” requires closed boundary shells."""
    assert part.solid.isClosed() and all(s.isClosed() for s in part.solid.Shells), (
        "The material boundary is open"
    )


@ALL_PARTS
def test_positive_volume(part):
    """“Positive-volume” requires a material solid after STEP exchange."""
    assert part.solid.Volume > 0, "The STEP solid has no positive material volume"


@ALL_PARTS
def test_original_stock_near_surrounding_faces(part):
    """“Preserve the surrounding shoulders, lead-ins and openings” on their material side."""
    check_occupancy(part, "normal_stock")


@ALL_PARTS
def test_original_void_near_surrounding_faces(part):
    """“Apply the rule in both material and void, including stock outside the affected bands.”"""
    check_occupancy(part, "normal_void")


@WITH_SEATS
def test_nominal_bores_stay_open(part):
    """“At signed … radial offsets −0.07” the original bore space remains open."""
    check_occupancy(part, "bore_void")


@WITH_SEATS
def test_added_diametral_clearance(part):
    """“Another 0.40 mm of diametral clearance” removes stock inside eligible radial bands."""
    check_occupancy(part, "allowance")


@WITH_SEATS
def test_material_behind_enlarged_walls(part):
    """“Strictly between r and r+0.20 mm” limits the depth of the wall relief."""
    check_occupancy(part, "backing")


@WITH_SEATS
def test_axial_engagement_at_trim_edges(part):
    """“Keep the same axial engagement” at axial offsets from the actual trim edges."""
    check_occupancy(part, "axial_trim")


@WITH_SEATS
def test_hoop_extent_at_trim_edges(part):
    """“Use the actual footprint, including disconnected pieces and internal boundaries.”"""
    check_occupancy(part, "hoop_trim")


@pytest.mark.parametrize("part", ["mount", "scope_base"], indirect=True, scope="module")
def test_angular_coverage_uses_union(part):
    """“The union of the azimuths … measures at least 180°” excludes partial walls."""
    check_occupancy(part, "partial_support")


@pytest.mark.parametrize("part", ["scope_base", "stage"], indirect=True, scope="module")
def test_radius_limits_preserve_excluded_bores(part):
    """“Nominal radius is in [1.25,10] mm” preserves fully covered walls outside the interval."""
    check_occupancy(part, "radius_excluded")


@ALL_PARTS
def test_single_delivered_cad_product(part):
    """“Write only your final deliverables” requires the revised STEP as the sole CAD product."""
    assert not part.extras, f"Additional CAD products: {part.extras}"


@ALL_PARTS
def test_exported_boundary_preserves_material_and_void(part):
    """“Apply the rule in both material and void, including stock outside the affected bands.”"""
    check_occupancy(part, "export_boundary")


@pytest.mark.parametrize("part", ["slide"], indirect=True, scope="module")
def test_nominal_part_without_relief(part):
    """“A part with no qualifying wall retains its nominal geometry” in its original frame."""
    check_occupancy(part, "nominal_material")


@pytest.mark.parametrize("part", ["stage"], indirect=True, scope="module")
def test_minimum_span_preserves_short_wall(part):
    """“The span from its lowest to highest occupied axial station is at least 1 mm.”"""
    check_occupancy(part, "span_excluded")
