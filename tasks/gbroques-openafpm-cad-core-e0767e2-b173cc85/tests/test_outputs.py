import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from functools import cache

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
GRADER = Path(__file__).parent / "grader"
sys.path.insert(0, str(GRADER))
PYTHON = "/opt/conda/bin/python"
PACKAGE = "openafpm_cad_core"
SOURCE = Path("/app/inputs") / PACKAGE
OUT = Path("/app/out")
# The omitted value is an independent API case even at the lower endpoint.
PRESETS = ("H Shape", "Star Shape", "T Shape")
# Cut lengths in mm. Omission retains the preset length.
CUTS = {
    "H Shape": {"0.10": 1080, "0.15": 2520, "0.14821": 1109, "default": None},
    "Star Shape": {"0.10": 3640, "0.15": 1560, "0.14273": 3511, "default": None},
    # At the long T cut the optional fraction is omitted (the lower endpoint).
    # The directly opened T assembly covers the default cut and fraction together.
    "T Shape": {"default": 1400, "0.15": 600},
}
LIVE = [(p, q) for p in PRESETS for q in CUTS[p]]
SAVED = ("T Shape", "saved")
EDITED = [("T Shape", "edit-0.1374"), ("T Shape", "edit-0.15")]
CONFIGURED = [*LIVE, *EDITED]
ALL = [SAVED, *CONFIGURED]
POSITION = 0.05
ANGLE = 0.01


def fraction(setting):
    value = setting[1].removeprefix("edit-")
    return 0.10 if value in {"saved", "default"} else float(value)


def assembly_entry(setting):
    if setting in EDITED:
        return "edited"
    return "saved" if setting == SAVED else "assembly"


def cases(settings):
    return pytest.mark.parametrize("setting", settings, scope="module", ids=[f"{p.split()[0]}-{q}" for p, q in settings])


@contextmanager
def work_directory():
    # The invoker can write its working directory but cannot replace that
    # directory or the independent reader through a writable parent.
    with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
        parent = Path(tmp)
        folder = parent / "work"
        folder.mkdir(mode=0o700)
        os.chown(folder, 65534, 65534)
        yield folder


def run_delivered(command, scratch, cwd=None):
    cwd = cwd or scratch
    env = {k: v for k, v in os.environ.items()
           if k in {"PATH", "LANG", "LC_ALL", "LD_LIBRARY_PATH", "QT_QPA_PLATFORM",
                    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "CONDA_PREFIX"}}
    env.update(HOME=str(cwd), TMPDIR=str(cwd), MPLCONFIGDIR=str(cwd))
    env["INVOKER_PARENT_PID"] = str(os.getpid())
    # Serialize access to this directory inode. Close the original input and
    # delivery paths to the child without changing any file in either tree.
    app = os.open("/app", os.O_DIRECTORY | os.O_RDONLY)
    try:
        fcntl.flock(app, fcntl.LOCK_EX)
        mode = stat.S_IMODE(os.fstat(app).st_mode)
        try:
            os.fchmod(app, mode & 0o700)
            scratch.parent.chmod(0o711)
            with (scratch.parent / "stdout").open("wb") as out, (scratch.parent / "stderr").open("wb") as err:
                return subprocess.run([sys.executable, "-I", "-B", str(GRADER / "invoke.py"), *command],
                                      cwd=cwd, env=env, stdin=subprocess.DEVNULL,
                                      stdout=out, stderr=err).returncode
        finally:
            # invoke.py has reaped this call's descendants. Close its private
            # tree before inspecting any outputs as root.
            try:
                scratch.parent.chmod(0o700)
                subprocess.run(["find", str(scratch), "!", "-type", "d", "!", "-type", "f", "-delete"], check=True)
            finally:
                os.fchmod(app, mode)
    finally:
        os.close(app)


def read_json(path):
    assert path.is_file(), f"Native measurement missing: {path.name}"
    with path.open() as file:
        return json.load(file)


def diagnostic(path):
    if not path.exists() or not stat.S_ISREG(path.lstat().st_mode):
        return ""
    with path.open("rb") as file:
        file.seek(max(0, path.stat().st_size - 3000))
        tail = file.read().decode(errors="replace").strip()
    # Lead with the exception so a one-line report identifies the failed call.
    # Keep its final frame and source line for the full pytest report.
    frame = tail.rfind('\n  File "')
    if frame >= 0:
        location, _, message = tail[frame:].strip().rpartition("\n")
        if location:
            return f"{message}\n{location}"
    return tail


def copy_for_call(source, scratch):
    target = scratch / "call"
    shutil.copytree(source, target, symlinks=True)
    for node in (target, *target.rglob("*")):
        os.chown(node, 65534, 65534, follow_symlinks=False)
    return target


def reader_modules(scratch):
    reader = scratch.parent / "reader"
    reader.mkdir(mode=0o755)
    for name in ("native.py", "frames.py", "solids.py", "native_fields.py"):
        shutil.copyfile(GRADER / name, reader / name)
    return reader


@pytest.fixture(scope="session")
def release(tmp_path_factory):
    source = OUT / PACKAGE
    assert source.is_dir() and not source.is_symlink(), f"Missing revised package at {source}"
    for p in source.rglob("*"):
        if p.is_symlink():
            assert p.resolve(strict=True).is_relative_to(source), f"Dependency escapes release: {p}"
    tree = tmp_path_factory.mktemp("release") / "delivery"
    shutil.copytree(OUT, tree, symlinks=True)
    return tree / PACKAGE


def native_reader(package, document, mode, folder, *args):
    # A new HOME, working directory and Python interpreter prevent the live
    # call from supplying measurements, reader modules or a FreeCAD user profile.
    with work_directory() as scratch:
        work = copy_for_call(package.parent, scratch)
        target = work / PACKAGE
        source = work / document
        assert source.is_relative_to(work) and ".." not in source.parts and source.is_file(), (
            f"Expected a native production assembly inside the release, received {document}")
        reader = reader_modules(scratch)
        record = scratch / "record"
        status = run_delivered([PYTHON, "-I", "-B", str(reader / "native.py"),
                                str(source), str(record), mode, str(target), *map(str, args)], scratch, cwd=work)
        if status:
            return None, f"Cannot inspect the native assembly ({status}): {diagnostic(scratch.parent/'stderr')}"
        assert record.is_dir(), "No assembly measurements returned from the native reader"
        shutil.copytree(record, folder, dirs_exist_ok=True)
    return ({} if mode in {"trial", "edited_trial"} else read_json(folder / "result.json")), ""


@pytest.fixture(scope="session")
def capture(release, tmp_path_factory, request):
    def run(setting, entry):
        preset, q = setting
        folder = tmp_path_factory.mktemp("captured")
        cut = CUTS[preset].get(q)
        if entry in {"saved", "edited"}:
            args = (fraction(setting),) if entry == "edited" else ()
            model, error = native_reader(release, PACKAGE + "/documents/WindTurbine.FCStd", entry, folder, *args)
            if model is not None:
                model["document"] = PACKAGE + "/documents/WindTurbine.FCStd"
            return (folder, model), error
        with work_directory() as scratch:
            work = copy_for_call(release.parent, scratch)
            package = work / PACKAGE
            # -c starts in the package parent as specified. The adapter imports
            # only stdlib before the public import sequence.
            reader = reader_modules(scratch)
            code = (GRADER / "capture.py").read_text() + "\ninvoke(*sys.argv[1:])"
            status = run_delivered([PYTHON, "-B", "-c", code, entry, preset, q.removeprefix("edit-"), str(cut or ""),
                                    str(reader), str(scratch)], scratch, cwd=work)
            if status:
                return None, (f"{entry} failed for {preset}, B={cut if cut is not None else 'default'}, "
                              f"q={q} ({status}): {diagnostic(scratch.parent/'stderr')}")
            data = read_json(scratch / "result.json")
            if entry == "tables":
                from comparator import illustration_inventory
                data["illustrations"] = illustration_inventory(data["tables"], package)
                return (folder, data), ""
            handoff = folder / "handoff"
            shutil.copytree(scratch / "mirror", handoff)
            shutil.copytree(scratch / "returned", folder / "returned")
        native, error = native_reader(handoff / PACKAGE, data["document"], "assembly", folder, data["master"])
        if error:
            return None, error
        data.update(native)
        return (folder, data), ""

    # The public calls remain serial and use fresh work trees. While pytest
    # measures one completed model, the worker can prepare the next model.
    jobs = {}
    with ThreadPoolExecutor(max_workers=1) as worker:
        def schedule(setting, entry):
            key = setting, entry
            if key not in jobs:
                jobs[key] = worker.submit(run, setting, entry)
            return jobs[key]

        for item in request.session.items:
            setting = getattr(getattr(item, "callspec", None), "params", {}).get("setting")
            if not setting:
                continue
            if "assembled" in item.fixturenames or "reference" in item.fixturenames:
                schedule(setting, assembly_entry(setting))
            if item.originalname in {"test_independent_fabrication_dimensions", "test_fabrication_illustrations"}:
                schedule(setting, "tables")

        def capture_result(setting, entry):
            result, error = schedule(setting, entry).result()
            assert not error and result is not None, error
            return result
        try:
            yield capture_result
        finally:
            worker.shutdown(wait=True, cancel_futures=True)


@pytest.fixture(scope="session")
def assembled(capture):
    return lambda setting: capture(setting, assembly_entry(setting))


@pytest.fixture(scope="session")
def reference(assembled, tmp_path_factory):
    membership_dir = tmp_path_factory.mktemp("membership")

    @cache
    def build(setting):
        preset, angle, width, cut = setting
        folder = tmp_path_factory.mktemp("native")
        shutil.copytree(SOURCE, folder / PACKAGE)
        env = dict(os.environ, HOME=str(folder), TMPDIR=str(folder), QT_QPA_PLATFORM="offscreen")
        with (folder / "stdout").open("wb") as out, (folder / "stderr").open("wb") as err:
            # This process contains only the supplied input and measurement code.
            # Cache a failed regeneration too, so dependent checks can report it.
            try:
                status = subprocess.run(
                    [PYTHON, "-P", "-B", str(GRADER / "reference.py"), preset,
                     str(angle), str(width), str(cut or ""), str(membership_dir / (preset + ".json"))],
                    cwd=folder, env=env, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                    timeout=300,
                ).returncode
            except subprocess.TimeoutExpired:
                err.write(b"\nThe independent comparator did not finish regeneration.\n")
                status = 1
        return folder, status

    @cache
    def at_released_settings(setting, before_rebuild=False):
        _, model = assembled(setting)
        preset = setting[0]
        angle = model["returned_angle"] if before_rebuild else model["angle"]
        width = model["returned_width"] if before_rebuild else model["width"]
        assert math.isfinite(angle) and math.isfinite(width), "The released stop settings must be finite"
        folder, status = build((preset, angle, width, CUTS[preset].get(setting[1])))
        assert status == 0, f"Cannot regenerate the comparator at {angle}, {width}: {diagnostic(folder/'stderr')}"
        return folder / "record", read_json(folder / "record" / "result.json")
    return at_released_settings


@pytest.fixture(scope="session")
def preserved(assembled, reference):
    from occurrences import compare_inventory
    @cache
    def result(setting):
        return compare_inventory(reference(setting)[0], assembled(setting)[0])
    return result


@pytest.fixture(scope="session")
def mechanism(assembled, reference, release):
    from motion import analyze, resolve_trial
    @cache
    def result(setting):
        folder, model = assembled(setting)
        original, data = reference(setting)
        selection = resolve_trial(original, folder)
        if selection is not None:
            if setting in EDITED:
                _, error = native_reader(release, model["document"], "edited_trial", folder,
                                         fraction(setting), *selection)
            else:
                source = release if setting == SAVED else folder / "handoff" / PACKAGE
                _, error = native_reader(source, model["document"], "trial", folder, *selection)
            assert not error, error
        radius = data["config"]["magnafpm"]["RotorDiameter"] / 2
        q = fraction(setting)
        metrics = analyze(original, folder, radius, q, model.get("motion"), boom_change=data["boom_change"])
        metrics["target"] = q * radius
        print(f"Measured {setting}: {metrics}; native limit {model['angle']} degrees")
        return metrics
    return result


@cases(ALL)
def test_settled_native_recompute(setting, assembled):
    """The native production graph recomputes to settled dimensions in each build."""
    _, model = assembled(setting)
    if setting in EDITED:
        assert not model["returned_invalid_objects"], f"Master edit did not settle: {model['returned_invalid_objects']}"
    assert not model["invalid_objects"], f"Native objects failed to recompute: {model['invalid_objects']}"


@cases(ALL)
def test_as_delivered_rest_assembly(setting, assembled, reference):
    """The live return, saved default and recomputed master edit each describe the fitted machine."""
    from occurrences import compare_inventory, material_deviation
    folder, model = assembled(setting)
    nominal = reference(setting, setting in EDITED)[0]
    observed = folder / ("opened" if setting == SAVED else "returned")
    result = compare_inventory(nominal, observed)
    assert result["count"] == result["expected_count"] and not result["invalid"], f"Incomplete rest assembly: {result}"
    assert result["angle_deg"] is not None and result["angle_deg"] <= ANGLE, f"Wrong rest datums: {result}"
    deviation, detail = material_deviation(nominal, observed)
    assert deviation <= POSITION, f"As-delivered material: {detail}"
    if setting in LIVE:
        assert abs(model["returned_angle"] - model["angle"]) <= ANGLE, "The unrounded limit changed during handoff"
        assert abs(model["returned_width"] - model["width"]) <= POSITION, "The fitted width changed during handoff"


@cases(ALL)
def test_as_delivered_furl_geometry(setting, assembled, reference):
    """The delivered or edited assembly meets the requested clearance before a forced rebuild."""
    from motion import analyze
    folder, model = assembled(setting)
    nominal, data = reference(setting, setting in EDITED)
    observed = folder / ("opened" if setting == SAVED else "returned")
    radius = data["config"]["magnafpm"]["RotorDiameter"] / 2
    q = fraction(setting)
    motion = model.get("returned_motion", model.get("motion"))
    m = analyze(nominal, observed, radius, q, motion, boom_change=data["boom_change"])
    limit = model.get("returned_angle", model["angle"])
    assert m["wall_contact"] and m["overlap"] > 0, f"Returned stop does not block on the yaw wall: {m}"
    assert abs(m["contact"] - limit) <= ANGLE, f"Incorrect as-delivered native limit: {m}"
    assert abs(m["clearance"] - q * radius) <= POSITION and m["minimum"] >= q * radius - POSITION, (
        f"Incorrect as-delivered vane clearance: {m}")
    assert m["motion_linear"] <= POSITION and m["motion_angular"] <= ANGLE, f"Wrong returned motion datums: {m}"


def test_single_released_turbine(release):
    """The final delivery contains the released turbine without a conflicting production copy."""
    from delivery import conflicting_assemblies
    root = Path("documents/WindTurbine.FCStd")
    extra = conflicting_assemblies(release.parent, release / root)
    assert not extra, f"Additional governed turbine assemblies in the delivery: {extra}"


@cases(ALL)
def test_moving_and_fixed_members(setting, mechanism):
    """The installed members retain the supplied moving and fixed relationships during furl."""
    assert mechanism(setting)["kinematic"], "Incorrect installed furling membership"


@cases(ALL)
def test_first_contact_on_outer_yaw_wall(setting, mechanism):
    """The finished high-end stop first contacts the outer cylindrical wall during the specified travel."""
    m = mechanism(setting)
    assert m["contact"] is not None and 0 <= m["contact"] <= 130, "No stop-to-pipe contact on the stroke"
    assert m["wall_contact"], f"First contact misses the native outer wall: {m}"


@cases(CONFIGURED)
def test_independent_fabrication_dimensions(setting, capture, reference):
    """A fresh table entry point reports the comparator's stock and dimensions at displayed precision."""
    from comparator import table_differences
    expected = reference(setting)[1]["tables"]
    _, data = capture(setting, "tables")
    delta = table_differences(expected, data["tables"])
    assert not delta, f"Fabrication table differences: {delta[:12]}"


@cases(CONFIGURED)
def test_fabrication_illustrations(setting, capture):
    """The complete relocated project supplies the drawings used by its fabrication tables."""
    _, data = capture(setting, "tables")
    assert not data["illustrations"], f"Unresolved fabrication illustrations: {data['illustrations']}"


@cases(ALL)
def test_interference_after_contact(setting, mechanism):
    """At 0.1 degrees beyond first contact the finished stop penetrates the finite yaw pipe."""
    m = mechanism(setting)
    assert m["overlap"] > 0, f"First contact does not block the stop: {m}"


@cases(ALL)
def test_finished_vane_clearance(setting, mechanism):
    """At contact the least vane clearance over the interchangeable boom lengths equals qR."""
    m = mechanism(setting)
    assert abs(m["clearance"] - m["target"]) <= POSITION, f"Incorrect full-furl clearance: {m}"


@cases(ALL)
def test_preceding_clearance(setting, mechanism):
    """Every service length retains the clearance reserve throughout the preceding stroke."""
    m = mechanism(setting)
    assert m["minimum"] >= m["target"] - POSITION, f"Insufficient clearance during the stroke: {m}"


@cases(ALL)
def test_native_travel(setting, mechanism, assembled):
    """The native unrounded limit agrees with the measured first blocking angle."""
    m = mechanism(setting)
    angle = assembled(setting)[1]["angle"]
    assert m["contact"] is not None and abs(angle - m["contact"]) <= ANGLE, f"{angle} vs {m['contact']} degrees"


@cases(CONFIGURED)
def test_published_travel(setting, assembled):
    """The motion export publishes the native travel rounded to two decimal places in degrees."""
    _, model = assembled(setting)
    maximum = model["motion"]["maximum_angle"]
    assert maximum == round(model["angle"], 2), "Incorrect rounded travel limit in settled model"
    if "returned_angle" in model:
        motion = model.get("returned_motion", model["motion"])
        assert motion["maximum_angle"] == round(model["returned_angle"], 2), "Incorrect rounded travel limit on return"


@cases(CONFIGURED)
def test_published_motion_datums(setting, mechanism):
    """The published transforms reproduce the installed parent, tail and outer-hinge datums."""
    m = mechanism(setting)
    assert m["motion_linear"] <= POSITION and m["motion_angular"] <= ANGLE, f"Motion datum mismatch: {m}"


@cases(ALL)
def test_installed_occurrences(setting, preserved):
    """Installed physical occurrences correspond one for one at the settled stop settings."""
    p = preserved(setting)
    assert p["count"] == p["expected_count"] and not p["missing_occurrences"], f"Installed occurrence differences: {p}"


@cases(ALL)
def test_valid_closed_parts(setting, preserved):
    """Every installed physical solid is valid and closed."""
    p = preserved(setting)
    assert p["count"] > 0 and not p["invalid"], f"Invalid or open installed solids: {p['invalid']}"


@cases(ALL)
def test_installed_material_boundaries(setting, preserved, reference, assembled):
    """Each installed material boundary agrees in both directions with the substituted comparator."""
    from occurrences import material_deviation
    p = preserved(setting)
    assert p["count"] == p["expected_count"] and not p["invalid"], f"Cannot pair closed material boundaries: {p}"
    deviation, detail = material_deviation(reference(setting)[0], assembled(setting)[0])
    assert deviation <= POSITION, f"Recomputed material: {detail}"


@cases(ALL)
def test_installed_datum_orientations(setting, preserved, mechanism):
    """Installed datum orientations and the datums that locate the moving mechanism are preserved."""
    p, m = preserved(setting), mechanism(setting)
    assert p["angle_deg"] is not None and p["angle_deg"] <= ANGLE, f"Installed datum orientation error: {p}"
    assert m["datum_angular"] <= ANGLE and m["datum_linear"] <= POSITION, f"Mechanism datum mismatch: {m}"


@cases(CONFIGURED)
def test_preset_inputs(setting, assembled, reference, capture):
    """The three named presets keep their original generator, furling and fabrication inputs."""
    _, result = capture(setting, "tables") if setting in EDITED else assembled(setting)
    groups = result["preset"]
    data = reference(setting)[1]
    base = dict(data["config"])
    base["furling"] = {**base["furling"], "BoomLength": data["nominal_boom"]}
    assert groups[1].get("VaneClearanceFraction", 0.10) == 0.10, "Changed default clearance fraction"
    for name, actual in zip(("magnafpm", "furling", "user"), groups, strict=True):
        for key, value in base[name].items():
            assert key in actual and actual[key] == value, (
                f"Preset input {name}.{key} should be {value!r}, received {actual.get(key)!r}")


@cases(ALL)
def test_relocated_native_links(setting, assembled):
    """The relocated production assembly resolves its linked native documents from its own copy."""
    data = assembled(setting)[1]
    files = data["external_documents"]
    assert not files, f"The assembly opened documents outside the relocated copy: {files}"

@cases(ALL)
def test_native_fabrication_and_preset(setting, assembled, reference):
    """The settled native build carries its preset and fitted settings through fabrication quantities."""
    from native_fields import native_differences
    q = fraction(setting)
    _, model = assembled(setting)
    differences = native_differences(model["native_data"], reference(setting)[1], SOURCE, q)
    if setting in EDITED:
        differences.extend(native_differences(model["returned_native_data"], reference(setting, True)[1], SOURCE, q))
    assert not differences, f"Native fabrication or preset differences: {differences[:12]}"
