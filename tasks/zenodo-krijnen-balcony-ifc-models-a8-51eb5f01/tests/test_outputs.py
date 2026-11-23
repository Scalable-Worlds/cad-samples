"""Checks for the cast-in-place concrete take-off: the delivered result for the input model and the
delivered program run on the input model and on other exports of the same project (held back with
their expected values)."""
import glob
import gzip
import json
import math
import os
import shutil
import stat
import subprocess
import tempfile
import time

import pytest

OUT = "/app/out"
SCRIPT = os.path.join(OUT, "pour_takeoff.py")
DELIVERED_RESULT = os.path.join(OUT, "a8_takeoff.json")
INPUT_MODEL = "/app/inputs/a8.ifc"
PYTHON = "/opt/conda/bin/python"

GRADER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "grader")
EXPECTED = os.path.join(GRADER, "expected.json")
HIDDEN_MODELS = {
    "relocated": os.path.join(GRADER, "models", "relocated.ifc.gz"),
    "regrouped": os.path.join(GRADER, "models", "regrouped.ifc.gz"),
    "reduced": os.path.join(GRADER, "models", "reduced.ifc.gz"),
    "rearranged": os.path.join(GRADER, "models", "rearranged.ifc.gz"),
    "combined": os.path.join(GRADER, "models", "combined.ifc.gz"),
    "restacked": os.path.join(GRADER, "models", "restacked.ifc.gz"),
    "resited": os.path.join(GRADER, "models", "resited.ifc.gz"),
}

KINDS = ("POER", "FUNDERINGSBALK", "VLOER", "OPSTORTING", "KOLOM", "KOLOMKOP", "WAND", "BALK", "DRUKLAAG")
QUANTITIES = ("volume_m3", "top_m2", "side_m2", "soffit_m2", "contact_m2")

UNPRIVILEGED = dict(user=65534, group=65534, extra_groups=[])


# ----------------------------------------------------------------------------------------------
# running the delivered program
# ----------------------------------------------------------------------------------------------
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
        except OSError:  # it ended meanwhile
            continue
        if "\nUid:\t65534\t" in status and "\nState:\tZ" not in status:
            return True
    return False


def run_delivered(command, scratch):
    """Run a delivered program as UID 65534 in `scratch`; returns its exit status."""
    env = {**os.environ, "HOME": scratch, "TMPDIR": scratch, "MPLCONFIGDIR": scratch}
    env.pop("PYTEST_CURRENT_TEST", None)  # the test id can contain an expected value
    with open(f"{scratch}/stdout", "wb") as out, open(f"{scratch}/stderr", "wb") as err:
        try:
            return subprocess.run(
                command, cwd=scratch, env=env, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                **UNPRIVILEGED,
            ).returncode
        finally:
            # Before anything the program wrote is read: `kill -9 -1`, run as UID 65534,
            # ends every process of that user; wait until they are gone.  Nothing the
            # program left behind is changed or deleted: what it wrote is read through
            # read_contained(), which never follows a link out of the program's own files.
            subprocess.run(["sh", "-c", "kill -9 -1"], env=env, **UNPRIVILEGED)
            while still_running():
                time.sleep(0.01)


def contained_file(path, roots):
    """The regular file that `path` names, directly or through symbolic links, when that file lies
    inside one of the directories `roots`; None when there is none or it lies anywhere else."""
    try:
        real = os.path.realpath(path, strict=True)
        info = os.lstat(real)
    except OSError:  # missing, dangling link, link loop
        return None
    if not stat.S_ISREG(info.st_mode):
        return None
    for root in roots:
        try:
            base = os.path.realpath(root, strict=True)
        except OSError:
            continue
        if real != base and os.path.commonpath([real, base]) == base:
            return real
    return None


def read_contained(path, roots, limit):
    """Bytes of the file `path` leads to, read with the rights of UID 65534, or None.

    The file must be a regular file of at most `limit` bytes inside one of `roots` (the program's
    scratch directory, the submission): a link to a file of the program's own is followed, a link that
    leaves those directories or that the program's user cannot read is refused."""
    real = contained_file(path, roots)
    if real is None:
        return None
    done = subprocess.run(
        ["head", "-c", str(limit + 1), "--", real],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, **UNPRIVILEGED,
    )
    if done.returncode != 0 or len(done.stdout) > limit:
        return None
    return done.stdout


def read_result(path, roots):
    """(parsed JSON, None) or (None, why the result cannot be read)."""
    raw = read_contained(path, roots, 4_000_000)
    if raw is None:
        return None, f"no readable result file at {path}"
    try:
        return json.loads(raw.decode("utf-8")), None
    except ValueError as exc:
        return None, f"the result at {path} is not JSON: {exc}"


class Outcome:
    """What one result file (delivered, or written by a run of the delivered program) holds."""

    def __init__(self, data=None, problem=None):
        self.data = data
        self.problem = problem


def run_on_model(name):
    """Runs the delivered program on the input model (at its own path) or on a held-back export."""
    if contained_file(SCRIPT, [OUT]) is None:
        return Outcome(problem=f"{SCRIPT} was not delivered")
    scratch = scratch_dir()
    try:
        if name == "a8":
            model = INPUT_MODEL
        else:
            model = os.path.join(scratch, "model.ifc")
            with gzip.open(HIDDEN_MODELS[name], "rb") as source, open(model, "wb") as target:
                shutil.copyfileobj(source, target)
            os.chmod(model, 0o644)
        result = os.path.join(scratch, "result.json")
        status = run_delivered([PYTHON, SCRIPT, model, result], scratch)
        if status != 0:
            tail = (read_contained(os.path.join(scratch, "stderr"), [scratch], 50_000_000) or b"")[-1500:]
            return Outcome(problem=f"the program exited with status {status}: {tail.decode('utf-8', 'replace')}")
        data, problem = read_result(result, [scratch, OUT])
        return Outcome(data=data, problem=problem)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


# ----------------------------------------------------------------------------------------------
# judging a result
# ----------------------------------------------------------------------------------------------
def number(data, kind, quantity):
    """The reported number, or None when the result does not hold one there."""
    try:
        value = data["kinds"][kind][quantity]
    except (KeyError, TypeError, IndexError):
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    return float(value)


def schema_problems(data):
    if not isinstance(data, dict) or not isinstance(data.get("kinds"), dict):
        return ['the result is not an object holding an object "kinds"']
    problems = []
    kinds = data["kinds"]
    if set(kinds) != set(KINDS):
        problems.append(f"kinds are {sorted(kinds)}, expected exactly {sorted(KINDS)}")
    for kind in KINDS:
        for quantity in QUANTITIES:
            if kind in kinds and number(data, kind, quantity) is None:
                problems.append(f"{kind}.{quantity} is not a finite number")
    return problems


def tolerance(quantity, exact):
    """Published accuracy: 0.5 % or 0.05 m3 for volumes, 1 % or 0.5 m2 for areas, whichever is larger."""
    if quantity == "volume_m3":
        return max(0.005 * abs(exact), 0.05)
    return max(0.01 * abs(exact), 0.5)


def exact_range(value):
    """(lowest, highest) exact value: one number, or a pair where the statement admits a range.

    The pair spans every reading the statement admits for the model, widened by the error bound of the
    computation the stored numbers come from: curved faces at their true shape (circles, not polygons)
    or faceted with a chordal error up to 1 mm; the slivers of negligible area whose layers lie between
    0.8 mm and 5 mm read as present or as absent; and, for a film of one kind thinner than 1 mm lying on
    concrete of another kind, the region owned by the kind behind the film (the film is a layer of
    concrete thinner than 1 mm, absent) or by the kind of the film (the concrete right behind the
    boundary of a body that is thick there).  The published tolerance applies at both ends, so a result
    within the tolerance of any of these values, the true-shape value among them, is accepted."""
    if isinstance(value, list):
        return float(min(value)), float(max(value))
    return float(value), float(value)


def assert_schema(outcome, what):
    assert outcome.problem is None, f"{what}: {outcome.problem}"
    problems = schema_problems(outcome.data)
    assert not problems, f"{what}: " + "; ".join(problems)


def assert_quantity(outcome, exact, quantity, what):
    assert outcome.problem is None, f"{what}: {outcome.problem}"
    misses = []
    for kind in KINDS:
        low, high = exact_range(exact[kind][quantity])
        got = number(outcome.data, kind, quantity)
        if got is None:
            misses.append(f"{kind}: no number reported")
        elif not low - tolerance(quantity, low) <= got <= high + tolerance(quantity, high):
            want = f"{low:.6f}" if low == high else f"{low:.6f} to {high:.6f}"
            allowed = tolerance(quantity, low if got < low else high)
            misses.append(f"{kind}: reported {got:.6f}, exact {want}, allowed deviation {allowed:.6f}")
    assert not misses, f"{what}, {quantity}: " + "; ".join(misses)


# ----------------------------------------------------------------------------------------------
# fixtures
# ----------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def expected():
    with open(EXPECTED) as file:
        return json.load(file)


@pytest.fixture(scope="module")
def delivered():
    data, problem = read_result(DELIVERED_RESULT, [OUT])
    return Outcome(data=data, problem=problem)


@pytest.fixture(scope="module")
def program_run():
    """Runs the delivered program once per model, on demand."""
    outcomes = {}

    def run(name):
        if name not in outcomes:
            outcomes[name] = run_on_model(name)
        return outcomes[name]

    return run


# ----------------------------------------------------------------------------------------------
# the delivered result for the input model
# ----------------------------------------------------------------------------------------------
def test_delivered_result_has_the_stated_structure(delivered):
    """Result: a8_takeoff.json is an object {"kinds": ...} with exactly the nine kinds, five numbers each."""
    assert_schema(delivered, "a8_takeoff.json")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_delivered_result_matches_the_input_model(delivered, expected, quantity):
    """Accuracy: every volume and area of a8_takeoff.json is within the published tolerance of a8.ifc's exact value."""
    assert_quantity(delivered, expected["a8"], quantity, "a8_takeoff.json")


# ----------------------------------------------------------------------------------------------
# the delivered program on the input model
# ----------------------------------------------------------------------------------------------
def test_input_model_run_writes_a_result(program_run):
    """Deliverables: run on /app/inputs/a8.ifc with a result path, the script exits with status 0 and writes the stated JSON."""
    assert_schema(program_run("a8"), "run on /app/inputs/a8.ifc")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_input_model_run_matches_the_input_model(program_run, expected, quantity):
    """Models/Accuracy: the script is run on /app/inputs/a8.ifc and its result is within the published tolerance."""
    assert_quantity(program_run("a8"), expected["a8"], quantity, "run on /app/inputs/a8.ifc")


# ----------------------------------------------------------------------------------------------
# the delivered program on other exports of the project
# ----------------------------------------------------------------------------------------------
def test_relocated_export_run_writes_a_result(program_run):
    """Deliverables: the script, given <model.ifc> <result.json>, exits with status 0 and writes the stated JSON."""
    assert_schema(program_run("relocated"), "export far from the origin, rotated, in metres")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_relocated_export_gives_the_same_takeoff(program_run, expected, quantity):
    """Models: within 1000 km of the origin, any rotation about the vertical axis, millimetres or metres."""
    assert_quantity(program_run("relocated"), expected["relocated"], quantity, "export far from the origin, rotated, in metres")


def test_regrouped_export_run_writes_a_result(program_run):
    """Deliverables: the script writes the stated JSON for an export whose parts are duplicated and subdivided."""
    assert_schema(program_run("regrouped"), "export with duplicated and subdivided parts")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_regrouped_export_gives_the_same_takeoff(program_run, expected, quantity):
    """Models: parts duplicated or subdivided; the union, its owners and its boundary partition do not change."""
    assert_quantity(program_run("regrouped"), expected["regrouped"], quantity, "export with duplicated and subdivided parts")


def test_reduced_export_run_writes_a_result(program_run):
    """Deliverables: the script writes the stated JSON for an export in which parts are missing."""
    assert_schema(program_run("reduced"), "export without some of the parts")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_reduced_export_takeoff(program_run, expected, quantity):
    """Models: exports differ in which parts exist; joints and contact turn into free surface where a neighbour is gone."""
    assert_quantity(program_run("reduced"), expected["reduced"], quantity, "export without some of the parts")


def test_rearranged_export_run_writes_a_result(program_run):
    """Deliverables: the script writes the stated JSON for an export in which parts are positioned differently."""
    assert_schema(program_run("rearranged"), "export with repositioned and added parts")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_rearranged_export_takeoff(program_run, expected, quantity):
    """Surface/Parts: layers under 1 mm are absent, nothing exists inside a void, the pour sequence owns overlaps."""
    assert_quantity(program_run("rearranged"), expected["rearranged"], quantity, "export with repositioned and added parts")


def test_combined_export_run_writes_a_result(program_run):
    """Deliverables: the script writes the stated JSON for an export that differs along all published axes."""
    assert_schema(program_run("combined"), "export that differs along all axes at once")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_combined_export_takeoff(program_run, expected, quantity):
    """Models: placement, rotation, unit, existing parts, positions, duplicates and subdivisions differ together."""
    assert_quantity(program_run("combined"), expected["combined"], quantity, "export that differs along all axes at once")


def test_restacked_export_run_writes_a_result(program_run):
    """Deliverables: the script writes the stated JSON for an export in which thin slices of parts are separate parts."""
    assert_schema(program_run("restacked"), "export with thin slices split off parts")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_restacked_export_takeoff(program_run, expected, quantity):
    """Surface: a stack of layers under 1 mm (a gap and a slice of concrete) is absent; what lies on either side touches."""
    assert_quantity(program_run("restacked"), expected["restacked"], quantity, "export with thin slices split off parts")


def test_resited_export_run_writes_a_result(program_run):
    """Deliverables: the script writes the stated JSON for a rotated export in metres with thin slices as separate parts."""
    assert_schema(program_run("resited"), "rotated export in metres with thin slices split off parts")


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_resited_export_takeoff(program_run, expected, quantity):
    """Surface/Models: stacks of layers under 1 mm are absent at any rotation, in metres, far from the origin."""
    assert_quantity(program_run("resited"), expected["resited"], quantity, "rotated export in metres with thin slices split off parts")
