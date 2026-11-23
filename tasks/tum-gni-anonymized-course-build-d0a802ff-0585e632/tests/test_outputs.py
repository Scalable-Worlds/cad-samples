"""Checks for the carrier operating chart of the roof trusses.

The delivered chart for the supplied structure is scored as it stands, and the
delivered script is run again on the supplied structure and on reissued
structures (the supplied model with its trusses replaced inside the published
reissue domain).  Every result is judged against a plane-frame reference that is
rebuilt from the IFC model in this process.
"""
import glob
import json
import math
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import time

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "grader"))

import frame_reference as reference  # noqa: E402
import knot_placement  # noqa: E402
import reissue_builder  # noqa: E402

PYTHON = "/opt/conda/bin/python"
OUT_DIR = "/app/out"
SCRIPT = os.path.join(OUT_DIR, "rate_trusses.py")
DELIVERED = os.path.join(OUT_DIR, "truss_ratings.json")
SUPPLIED = "/app/inputs/model_7_structure.ifc"
CASES_FILE = os.path.join(HERE, "grader", "reissue_cases.json")

REISSUES = ["short", "long", "tie", "joints", "sections", "offsets"]
RUNS = ["supplied"] + REISSUES
SOURCES = ["delivered"] + RUNS

# published tolerances
POSITION_TOL = 1.0          # mm, world positions and stations
AT_TOLERANCE = 1e-6         # mm: round-off allowed on a position or station that sits at its tolerance
AREA_TOL = 0.01
INERTIA_TOL = 0.02
FIBRE_TOL = 0.25            # mm
BAND_LOW, BAND_HIGH = 0.99, 1.001
WITNESS_TOL = 0.015
ROUND_OFF = 1e-9
# R(s) itself is only defined to about this relative accuracy: joints coincide within
# 0.05 mm and exact frame analyses that merge them differently differ by this much.
MODEL_NOISE = 2e-5

# How R(s) is followed along the chord.  These are properties of the reference, not limits on a
# curve: every knot and every stretch between two knots is judged, whatever their number or spacing.
GRID_STEP = 0.5             # mm: R(s) is evaluated at least this densely
FINEST = 2.0 ** -10         # mm: a stretch this short is not divided further
SHORTEST = 1e-6             # mm: the slope of R(s) is not taken over a shorter stretch
TURN_SAFETY = 4.0           # margin on the estimated departure of R(s) from a chord
SETTLED = 0.25 * MODEL_NOISE  # R(s) is not followed more closely than this share of itself

UNPRIVILEGED = dict(user=65534, group=65534, extra_groups=[])
PLAIN_ENV = {"PATH": "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}


# --------------------------------------------------------------------------
# running the delivered script
# --------------------------------------------------------------------------
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


def read_as_submission(path):
    """Bytes of a file as UID 65534 reads it, or None when it is not a regular file that this user can read.

    Symbolic links are followed, but only as far as the submission itself may read: a link
    into the tests or their logs gives nothing.  Called only while no delivered program is
    running."""
    try:
        done = subprocess.run(["sh", "-c", 'test -f "$1" && exec cat -- "$1"', "sh", path], cwd="/", env=PLAIN_ENV,
                              stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                              **UNPRIVILEGED)
    except OSError:
        return None
    return done.stdout if done.returncode == 0 else None


def directory_for_submission(path):
    """Whether UID 65534 finds a directory at `path` (links followed as that user)."""
    try:
        done = subprocess.run(["sh", "-c", 'test -d "$1"', "sh", path], cwd="/", env=PLAIN_ENV,
                              stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                              **UNPRIVILEGED)
    except OSError:
        return False
    return done.returncode == 0


def run_delivered(command, scratch, collect):
    """Run a delivered program as UID 65534 in `scratch`; returns its exit status and collect().

    `collect` takes what the program wrote.  It is called once the program and everything it
    started has ended, and it reads as UID 65534 (read_as_submission), so the result may be
    reached through symbolic links of the program's own making.  Only then are the symbolic
    links of the scratch directory deleted, before this process opens anything there itself.

    Nothing under /app/out is altered here, so every run meets the submission as it was
    delivered, and nothing below opens a path under /app/out once a delivered program has
    run."""
    env = {**os.environ, "HOME": scratch, "TMPDIR": scratch, "MPLCONFIGDIR": scratch}
    env.pop("PYTEST_CURRENT_TEST", None)  # the test id can contain an expected value
    collected = None
    with open(f"{scratch}/stdout", "wb") as out, open(f"{scratch}/stderr", "wb") as err:
        try:
            status = subprocess.run(
                command, cwd=scratch, env=env, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                **UNPRIVILEGED,
            ).returncode
        finally:
            # `kill -9 -1`, run as UID 65534, ends every process of that user; wait until
            # they are gone, so that nothing can change the files any more.
            subprocess.run(["sh", "-c", "kill -9 -1"], env=env, **UNPRIVILEGED)
            while still_running():
                time.sleep(0.01)
            try:
                collected = collect()
            finally:
                # through a link this process would read a protected file
                subprocess.run(["find", scratch, "-type", "l", "-delete"], check=True)
    return status, collected


def read_submitted(path):
    """Bytes of a delivered file, or None when there is none.

    Only called before any delivered program has run.  A regular file is opened
    without following links; a link is read as UID 65534, which yields only what the
    submission itself is allowed to read."""
    try:
        if stat.S_ISREG(os.lstat(path).st_mode):
            with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW), "rb") as file:
                return file.read()
    except OSError:
        return None
    return read_as_submission(path)


def _tail(path, limit=1500):
    try:
        with open(path, "rb") as file:
            data = file.read()
    except OSError:
        return ""
    return data[-limit:].decode("utf-8", "replace")


class Workspace(object):
    """Lazily built inputs, reference models and results, shared by all tests."""

    def __init__(self):
        # The delivered chart is read once, before any delivered program runs.
        self.delivered_bytes = read_submitted(DELIVERED)
        with open(CASES_FILE) as file:
            self.cases = json.load(file)
        self.vault = tempfile.mkdtemp(dir="/tmp")      # root only: master copies of the reissues
        self._models = {}
        self._references = {}
        self._runs = {}
        self._documents = {}
        self._ratings = {}
        self._followed = {}

    # ----------------------------------------------------------- inputs
    def model_path(self, source):
        if source in ("delivered", "supplied"):
            return SUPPLIED
        if source not in self._models:
            assert source in self.cases, "unknown reissue %s" % source
            path = os.path.join(self.vault, "%s.ifc" % source)
            reissue_builder.build_reissue(SUPPLIED, self.cases[source], path)
            self._models[source] = path
        return self._models[source]

    def trusses(self, source):
        path = self.model_path(source)
        if path not in self._references:
            self._references[path] = reference.load_trusses(path)
        return self._references[path]

    # ------------------------------------------------------------- runs
    def run(self, name):
        """Run the delivered script once on the named structure."""
        if name in self._runs:
            return self._runs[name]
        result = {"ran": False, "status": None, "out_dir_exists": False, "bytes": None, "note": ""}
        self._runs[name] = result
        if not os.path.isfile(SCRIPT):
            result["note"] = "%s is missing" % SCRIPT
            return result
        scratch = scratch_dir()
        if name == "supplied":
            structure = SUPPLIED
        else:
            structure = os.path.join(scratch, "structure_%s.ifc" % name)
            shutil.copyfile(self.model_path(name), structure)
            os.chmod(structure, 0o644)
        out_dir = os.path.join(scratch, "ratings")
        target = os.path.join(out_dir, "truss_ratings.json")
        result["status"], (result["out_dir_exists"], result["bytes"]) = run_delivered(
            [PYTHON, SCRIPT, "--structure", structure, "--out", out_dir], scratch,
            lambda: (directory_for_submission(out_dir), read_as_submission(target)))
        result["ran"] = True
        result["note"] = "stderr tail: %s" % _tail(os.path.join(scratch, "stderr"))
        return result

    def document(self, source):
        """Parsed chart of a source; fails the calling test when it is unavailable."""
        if source not in self._documents:
            if source == "delivered":
                assert self.delivered_bytes is not None, "%s is missing" % DELIVERED
                raw = self.delivered_bytes
            else:
                result = self.run(source)
                assert result["ran"], result["note"]
                assert result["status"] == 0, "script exited with status %s on %s; %s" % (
                    result["status"], source, result["note"])
                assert result["bytes"] is not None, "no truss_ratings.json in the --out directory for %s" % source
                raw = result["bytes"]
            try:
                parsed = json.loads(raw.decode("utf-8"))
            except (ValueError, UnicodeDecodeError) as exc:
                parsed = exc
            self._documents[source] = parsed
        parsed = self._documents[source]
        assert not isinstance(parsed, Exception), "truss_ratings.json of %s is not valid JSON: %s" % (source, parsed)
        return parsed

    def pairs(self, source):
        """(reference truss, delivered record) for every truss of the structure."""
        doc = self.document(source)
        assert isinstance(doc, dict) and isinstance(doc.get("trusses"), list), "no \"trusses\" array"
        records = {}
        for record in doc["trusses"]:
            assert isinstance(record, dict), "a truss entry is not an object"
            records.setdefault(record.get("guid"), []).append(record)
        out = []
        for truss in self.trusses(source):
            found = records.get(truss.guid, [])
            assert len(found) == 1, "assembly %s has %d records" % (truss.guid, len(found))
            out.append((truss, found[0]))
        return out

    def rating(self, source, truss, stations):
        """Reference R(s) in kN at the given stations; a station is evaluated once per truss."""
        stations = np.asarray(stations, dtype=float)
        key = (self.model_path(source), truss.guid)
        known, values = self._ratings.get(key, (np.empty(0), np.empty(0)))
        fresh = np.unique(stations)
        if known.size:
            fresh = fresh[known[np.minimum(np.searchsorted(known, fresh), known.size - 1)] != fresh]
        if fresh.size:
            known = np.concatenate([known, fresh])
            values = np.concatenate([values, truss.rating(fresh) / 1000.0])
            order = np.argsort(known, kind="stable")
            known, values = known[order], values[order]
            self._ratings[key] = (known, values)
        return values[np.searchsorted(known, stations)]

    def followed(self, source, truss):
        """R(s) of a truss as followed along 0..L-1100 (followed_rating), worked out once per truss."""
        key = (self.model_path(source), truss.guid)
        if key not in self._followed:
            self._followed[key] = followed_rating(truss, lambda stations: self.rating(source, truss, stations))
        return self._followed[key]


@pytest.fixture(scope="session")
def ws():
    return Workspace()


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def is_triple(value):
    return isinstance(value, list) and len(value) == 3 and all(is_number(v) for v in value)


def members_by_guid(truss, record):
    members = record.get("members")
    assert isinstance(members, list), "truss %s: \"members\" is not an array" % truss.guid
    table = {}
    for member in members:
        assert isinstance(member, dict), "truss %s: a member record is not an object" % truss.guid
        table[member.get("guid")] = member
    return table


def curve_arrays(truss, record):
    """Validated knots of a curve: stations (mm) and loads (kN)."""
    curve = record.get("curve")
    assert isinstance(curve, list) and len(curve) >= 2, "truss %s: curve needs at least two knots" % truss.guid
    for knot in curve:
        assert isinstance(knot, list) and len(knot) == 2 and is_number(knot[0]) and is_number(knot[1]), (
            "truss %s: curve knot %r is not a [s_mm, load_kN] pair" % (truss.guid, knot))
    stations = np.array([knot[0] for knot in curve], dtype=float)
    loads = np.array([knot[1] for knot in curve], dtype=float)
    assert np.all(np.diff(stations) > 0.0), "truss %s: curve stations are not strictly increasing" % truss.guid
    return stations, loads


def followed_rating(truss, rate):
    """Stations of 0..L-1100 and R(s) there, close enough for straight lines between them to follow R(s).

    R(s) is evaluated on a grid of GRID_STEP and wherever a wheel or an edge of
    the unit exclusion zones crosses a joint or support. Changes of governing response
    can occur between those events. Track them using departure estimates from the
    changes of slope at both ends of each evaluated stretch.  Every stretch whose bound exceeds
    SETTLED is halved, down to FINEST.

    Returns the stations, R(s) there in kN and the bound of every stretch between them."""
    s_max = truss.s_max
    events = np.concatenate([truss.xb - wheel - edge for wheel in (0., reference.WHEELBASE)
                             for edge in (0., -reference.CLEARANCE, reference.CLEARANCE)])
    grid = GRID_STEP * np.arange(0, math.floor(s_max / GRID_STEP) + 1)
    stations = np.concatenate([events, grid, [0.0, s_max]])
    stations = np.unique(stations[(stations >= 0.0) & (stations <= s_max)])
    while True:
        rated = rate(stations)
        width = np.diff(stations)
        depart = np.zeros(width.size)
        left = np.nonzero(width > SHORTEST)[0]        # long enough for a slope free of round-off
        if left.size >= 2:
            turn = np.abs(np.diff((rated[left + 1] - rated[left]) / width[left]))
            depart[left] = TURN_SAFETY * 0.25 * width[left] * (
                np.concatenate([turn[:1], turn]) + np.concatenate([turn, turn[-1:]]))
        halve = np.nonzero((width > FINEST) & (depart > SETTLED * np.maximum(rated[:-1], rated[1:])))[0]
        if halve.size == 0:
            return stations, rated, depart
        stations = np.union1d(stations, 0.5 * (stations[halve] + stations[halve + 1]))


def band_fault(truss, knots, loads, followed):
    """Why a curve is not inside the band 0.99 R .. 1.001 R, as a message, or None when it is.

    Stations carry 1 mm, so the knots of a compliant curve may stand up to POSITION_TOL
    from their analytical stations, each by its own amount.  A curve is inside the band
    when its knots can be given stations t such that

      * every t is within POSITION_TOL of the reported station of its knot,
      * the knots keep their order,
      * the first knot stands at or before 0 and the last one at or beyond L - 1100, and
      * the curve that is linear between the knots at these stations has its ordinate
        between 0.99 R(s) and 1.001 R(s) at every station s of 0..L-1100.

    Whether such stations exist is decided for all knots and all stretches between them
    together, exactly and without sampling the curve (knot_placement).  R(s) enters as
    `followed` (see followed_rating): straight between the evaluated stations, with the
    band widened by MODEL_NOISE and by the bound on the departure of R(s) from these
    straight lines.  The stations that were found are checked once more, directly, before
    the curve is accepted."""
    s_max = truss.s_max
    if not max(0.0, float(knots[0])) < min(s_max, float(knots[-1])):
        return "the curve has no extent on 0..%.3f mm" % s_max
    reach = POSITION_TOL + AT_TOLERANCE
    stations, rated, depart = followed
    doubt = np.maximum(np.concatenate([[0.0], depart]), np.concatenate([depart, [0.0]]))
    # The inner optimisation returns an attained upper capacity. Its certified lower
    # capacity is available independently of the joint-merging and interpolation allowances.
    from suspended_load import RELATIVE_GAP
    lower = BAND_LOW * (1.0 - MODEL_NOISE) * np.maximum((1.0 - RELATIVE_GAP) * rated - doubt, 0.0)
    upper = BAND_HIGH * (1.0 + MODEL_NOISE) * (rated + doubt)
    stuck, placed = knot_placement.place(knots, loads, stations, lower, upper, reach)
    if stuck < 0:
        beyond = knot_placement.excess(placed, knots, loads, stations, lower, upper, reach)
        if beyond > ROUND_OFF:
            return "the stations found for the knots do not hold when checked directly (off by %.3g)" % beyond
        return None
    if stuck == 0:
        return "the first knot, reported at s = %.4f mm, cannot stand at or before station 0" % knots[0]
    # what the curve does as reported, around the knot that cannot be placed
    a, b = max(0.0, float(knots[stuck - 1])), min(s_max, float(knots[stuck]))
    a, b = min(a, s_max), max(b, 0.0)
    along = np.unique(np.concatenate([stations[(stations >= a) & (stations <= b)], [a, b]]))
    charted = np.interp(along, knots, loads)
    capacity = np.interp(along, stations, rated)
    share = np.divide(charted, capacity, out=np.full(along.size, np.inf), where=capacity > 0.0)
    share[(capacity <= 0.0) & (charted == 0.0)] = 1.0
    worst = int(np.argmax(np.maximum(share / BAND_HIGH, BAND_LOW / np.maximum(share, 1e-300))))
    ends = np.clip([knots[stuck] - POSITION_TOL, knots[stuck] + POSITION_TOL], 0.0, s_max)
    near = np.concatenate([rated[(stations >= ends[0]) & (stations <= ends[1])], np.interp(ends, stations, rated)])
    return ("no stations within 1 mm of the reported ones, in the same order, keep the curve inside the band: "
            "knots 0..%d can be placed, knot %d (reported at s = %.4f mm with %.4f kN) cannot follow them.  "
            "R(s) within 1 mm of that knot is %.4f..%.4f kN.  As reported, the curve is at %.5f..%.5f of R(s) "
            "from s = %.4f to %.4f mm, worst at s = %.4f mm: charted %.4f kN, R(s) = %.4f kN" % (
                stuck - 1, stuck, knots[stuck], loads[stuck], near.min(), near.max(), share.min(), share.max(),
                a, b, along[worst], charted[worst], capacity[worst]))


def distance_to_segment(point, a, b):
    ab = b - a
    t = min(max(float((point - a) @ ab) / float(ab @ ab), 0.0), 1.0)
    return float(np.linalg.norm(point - (a + t * ab)))


def witness_of(truss, record):
    witness = record.get("witness")
    assert isinstance(witness, dict), "truss %s: witness is not an object" % truss.guid
    assert is_number(witness.get("s_mm")), "truss %s: witness s_mm is not a number" % truss.guid
    assert witness.get("heavy_at") in ("LOW", "HIGH"), "truss %s: witness heavy_at %r" % (
        truss.guid, witness.get("heavy_at"))
    assert witness.get("mode") in reference.MODES, "truss %s: witness mode %r" % (truss.guid, witness.get("mode"))
    assert isinstance(witness.get("member"), str), "truss %s: witness member is not a GlobalId" % truss.guid
    assert is_triple(witness.get("xyz_mm")), "truss %s: witness xyz_mm is not three numbers" % truss.guid
    assert -POSITION_TOL - AT_TOLERANCE <= witness["s_mm"] <= truss.s_max + POSITION_TOL + AT_TOLERANCE, (
        "truss %s: witness station %.3f is outside 0..%.3f" % (truss.guid, witness["s_mm"], truss.s_max))
    return witness


def witness_member(truss, witness):
    """The reference member named by a witness."""
    for member in truss.members:
        if member["guid"] == witness["member"]:
            return member
    raise AssertionError("truss %s: witness member %r is not one of its members" % (truss.guid, witness["member"]))


# --------------------------------------------------------------------------
# the script
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", RUNS)
def test_script_writes_ratings(ws, name):
    """The stated command works from another working directory, creates --out and writes truss_ratings.json there."""
    result = ws.run(name)
    assert result["ran"], result["note"]
    assert result["status"] == 0, "exit status %s; %s" % (result["status"], result["note"])
    assert result["out_dir_exists"], "the --out directory was not created"
    assert result["bytes"] is not None, "truss_ratings.json was not written directly in --out"
    assert isinstance(ws.document(name), dict), "truss_ratings.json is not a JSON object"


# --------------------------------------------------------------------------
# structure of the report
# --------------------------------------------------------------------------
@pytest.mark.parametrize("source", SOURCES)
def test_report_schema(ws, source):
    """The JSON object has the one key "trusses" and each entry carries guid, members, curve and witness as specified."""
    doc = ws.document(source)
    assert isinstance(doc, dict), "top level is not an object"
    assert set(doc) == {"trusses"}, "top-level keys are %s" % sorted(doc)
    assert isinstance(doc["trusses"], list) and doc["trusses"], "\"trusses\" is not a non-empty array"
    for record in doc["trusses"]:
        assert isinstance(record, dict), "a truss entry is not an object"
        for key in ("guid", "members", "curve", "witness"):
            assert key in record, "truss entry without %r" % key
        assert isinstance(record["guid"], str), "truss guid is not a string"
        assert isinstance(record["members"], list) and record["members"], "members is not a non-empty array"
        for member in record["members"]:
            assert isinstance(member, dict), "member record is not an object"
            for key in ("guid", "ends_mm", "A_mm2", "I_mm4", "c_mm"):
                assert key in member, "member record without %r" % key
            assert isinstance(member["guid"], str), "member guid is not a string"
            ends = member["ends_mm"]
            assert isinstance(ends, list) and len(ends) == 2 and all(is_triple(e) for e in ends), (
                "ends_mm of %s is not two XYZ triples" % member["guid"])
            assert is_number(member["A_mm2"]) and member["A_mm2"] > 0, "A_mm2 of %s" % member["guid"]
            assert is_number(member["I_mm4"]) and member["I_mm4"] > 0, "I_mm4 of %s" % member["guid"]
            assert member["c_mm"] is None or is_number(member["c_mm"]), "c_mm of %s" % member["guid"]
        curve = record["curve"]
        assert isinstance(curve, list) and len(curve) >= 2, "curve of %s needs at least two knots" % record["guid"]
        for knot in curve:
            assert isinstance(knot, list) and len(knot) == 2 and all(is_number(v) for v in knot), (
                "curve knot %r of %s is not a [s_mm, load_kN] pair" % (knot, record["guid"]))
        witness = record["witness"]
        assert isinstance(witness, dict), "witness of %s is not an object" % record["guid"]
        for key in ("s_mm", "unit_s_mm", "heavy_at", "mode", "member", "xyz_mm"):
            assert key in witness, "witness of %s without %r" % (record["guid"], key)
        assert is_number(witness["s_mm"]), "witness s_mm of %s" % record["guid"]
        assert is_number(witness["unit_s_mm"]), "witness unit_s_mm of %s" % record["guid"]
        assert witness["heavy_at"] in ("LOW", "HIGH"), "witness heavy_at of %s" % record["guid"]
        assert witness["mode"] in reference.MODES, "witness mode of %s" % record["guid"]
        assert isinstance(witness["member"], str), "witness member of %s" % record["guid"]
        assert is_triple(witness["xyz_mm"]), "witness xyz_mm of %s" % record["guid"]


@pytest.mark.parametrize("source", SOURCES)
def test_one_record_per_assembly(ws, source):
    """"trusses" holds exactly one object per TRUSS assembly, keyed by the assembly GlobalId."""
    doc = ws.document(source)
    assert isinstance(doc, dict) and isinstance(doc.get("trusses"), list), "no \"trusses\" array"
    got = [record.get("guid") if isinstance(record, dict) else None for record in doc["trusses"]]
    expected = sorted(truss.guid for truss in ws.trusses(source))
    assert sorted(str(g) for g in got) == expected, "assembly GlobalIds differ: got %s" % sorted(str(g) for g in got)


@pytest.mark.parametrize("source", SOURCES)
def test_one_record_per_member(ws, source):
    """"members" holds exactly one record per aggregated IfcMember, keyed by the member GlobalId."""
    for truss, record in ws.pairs(source):
        members = record.get("members")
        assert isinstance(members, list), "truss %s: members is not an array" % truss.guid
        got = sorted(str(m.get("guid")) if isinstance(m, dict) else "?" for m in members)
        expected = sorted(m["guid"] for m in truss.members)
        assert got == expected, "truss %s: member GlobalIds differ (%d reported, %d aggregated)" % (
            truss.guid, len(got), len(expected))


# --------------------------------------------------------------------------
# member geometry and sections
# --------------------------------------------------------------------------
@pytest.mark.parametrize("source", SOURCES)
def test_member_axis_ends(ws, source):
    """ends_mm gives the two project-world ends of each placed Axis segment within 1 mm, in either order."""
    for truss, record in ws.pairs(source):
        table = members_by_guid(truss, record)
        for member in truss.members:
            got = table.get(member["guid"])
            assert got is not None, "truss %s: member %s is missing" % (truss.guid, member["guid"])
            ends = got.get("ends_mm")
            assert isinstance(ends, list) and len(ends) == 2 and all(is_triple(e) for e in ends), (
                "member %s: ends_mm is not two XYZ triples" % member["guid"])
            a, b = np.array(ends[0], dtype=float), np.array(ends[1], dtype=float)
            direct = max(np.linalg.norm(a - member["p1"]), np.linalg.norm(b - member["p2"]))
            swapped = max(np.linalg.norm(a - member["p2"]), np.linalg.norm(b - member["p1"]))
            assert min(direct, swapped) <= POSITION_TOL + AT_TOLERANCE, (
                "member %s: axis ends are %.3f mm off" % (member["guid"], min(direct, swapped)))


@pytest.mark.parametrize("source", SOURCES)
def test_section_area(ws, source):
    """A_mm2 is the steel area of the placed Body section at mid-length, within 1%."""
    for truss, record in ws.pairs(source):
        table = members_by_guid(truss, record)
        for member in truss.members:
            got = table.get(member["guid"])
            assert got is not None, "truss %s: member %s is missing" % (truss.guid, member["guid"])
            assert is_number(got.get("A_mm2")), "member %s: A_mm2 is not a number" % member["guid"]
            assert abs(got["A_mm2"] / member["A"] - 1.0) <= AREA_TOL + ROUND_OFF, (
                "member %s: A_mm2 %.3f is off by %.2f%%" % (
                    member["guid"], got["A_mm2"], 100.0 * (got["A_mm2"] / member["A"] - 1.0)))


@pytest.mark.parametrize("source", SOURCES)
def test_section_second_moment(ws, source):
    """I_mm4 is I about the truss-plane normal for a chord and the smaller principal I for a web, within 2%."""
    for truss, record in ws.pairs(source):
        table = members_by_guid(truss, record)
        for member in truss.members:
            got = table.get(member["guid"])
            assert got is not None, "truss %s: member %s is missing" % (truss.guid, member["guid"])
            assert is_number(got.get("I_mm4")), "member %s: I_mm4 is not a number" % member["guid"]
            assert abs(got["I_mm4"] / member["I"] - 1.0) <= INERTIA_TOL + ROUND_OFF, (
                "%s %s: I_mm4 %.6g is off by %.2f%%" % (
                    member["kind"], member["guid"], got["I_mm4"], 100.0 * (got["I_mm4"] / member["I"] - 1.0)))


@pytest.mark.parametrize("source", SOURCES)
def test_section_fibre_distance(ws, source):
    """c_mm is the distance to the most distant in-plane fibre for a chord (0.25 mm) and null for a web."""
    for truss, record in ws.pairs(source):
        table = members_by_guid(truss, record)
        for member in truss.members:
            got = table.get(member["guid"])
            assert got is not None, "truss %s: member %s is missing" % (truss.guid, member["guid"])
            assert "c_mm" in got, "member %s: no c_mm" % member["guid"]
            if member["kind"] == "web":
                assert got["c_mm"] is None, "web %s: c_mm must be null" % member["guid"]
            else:
                assert is_number(got["c_mm"]), "chord %s: c_mm is not a number" % member["guid"]
                assert abs(got["c_mm"] - member["c"]) <= FIBRE_TOL + ROUND_OFF, (
                    "chord %s: c_mm %.3f is %.3f mm off" % (member["guid"], got["c_mm"], got["c_mm"] - member["c"]))


# --------------------------------------------------------------------------
# the chart
# --------------------------------------------------------------------------
@pytest.mark.parametrize("source", SOURCES)
def test_curve_knots_cover_the_chord(ws, source):
    """The curve has strictly increasing knots of finite [s_mm, load_kN] covering stations 0 through L - 1100 (1 mm)."""
    for truss, record in ws.pairs(source):
        stations, _ = curve_arrays(truss, record)
        assert stations[0] <= POSITION_TOL + AT_TOLERANCE, "truss %s: first knot at %.3f does not cover station 0" % (
            truss.guid, stations[0])
        assert stations[-1] >= truss.s_max - POSITION_TOL - AT_TOLERANCE, (
            "truss %s: last knot at %.3f does not reach %.3f" % (truss.guid, stations[-1], truss.s_max))


@pytest.mark.parametrize("source", SOURCES)
def test_curve_within_band(ws, source):
    """At every station the chart lies within 0.99..1.001 of the capacity for both orientations and all legal unit positions."""
    for truss, record in ws.pairs(source):
        knots, loads = curve_arrays(truss, record)
        assert knots[0] <= POSITION_TOL + AT_TOLERANCE and knots[-1] >= truss.s_max - POSITION_TOL - AT_TOLERANCE, (
            "truss %s: curve does not cover 0..%.3f" % (truss.guid, truss.s_max))
        fault = band_fault(truss, knots, loads, ws.followed(source, truss))
        assert fault is None, "truss %s: %s" % (truss.guid, fault)


# --------------------------------------------------------------------------
# the witness
# --------------------------------------------------------------------------
@pytest.mark.parametrize("source", SOURCES)
def test_witness_at_curve_minimum(ws, source):
    """The witness station is within 1 mm of a minimum of the delivered curve."""
    for truss, record in ws.pairs(source):
        knots, loads = curve_arrays(truss, record)
        witness = witness_of(truss, record)
        lo, hi = max(0.0, knots[0]), min(truss.s_max, knots[-1])
        assert lo < hi, "truss %s: curve has no extent on the chord" % truss.guid
        inside = knots[(knots >= lo) & (knots <= hi)]
        domain = np.unique(np.concatenate([inside, [lo, hi]]))
        minimum = float(np.min(np.interp(domain, knots, loads)))
        # The curve is linear between its knots, so its lowest ordinate within 1 mm of the
        # witness is at a knot or at an end of that range.  Minima are equal up to round-off.
        near = np.concatenate([[witness["s_mm"] - POSITION_TOL, witness["s_mm"], witness["s_mm"] + POSITION_TOL],
                               inside[np.abs(inside - witness["s_mm"]) <= POSITION_TOL + AT_TOLERANCE]])
        at_witness = float(np.min(np.interp(np.clip(near, lo, hi), knots, loads)))
        assert at_witness <= minimum + ROUND_OFF * max(abs(minimum), 1e-3), (
            "truss %s: within 1 mm of the witness station %.3f mm the curve does not go below %.9g kN, "
            "but its minimum is %.9g kN" % (truss.guid, witness["s_mm"], at_witness, minimum))


@pytest.mark.parametrize("source", SOURCES)
def test_witness_names_member_and_point(ws, source):
    """The witness names a member of the truss that suits its mode, with the chord axis point or the web midpoint (1 mm)."""
    for truss, record in ws.pairs(source):
        witness = witness_of(truss, record)
        member = witness_member(truss, witness)
        point = np.array(witness["xyz_mm"], dtype=float)
        if witness["mode"] == "SAG":
            assert member is truss.low, "truss %s: SAG must name the lower chord" % truss.guid
        if witness["mode"] == "WEB_BUCKLING":
            assert member["kind"] == "web", "truss %s: WEB_BUCKLING must name a web" % truss.guid
        if member["kind"] == "web":
            off = float(np.linalg.norm(point - 0.5 * (member["p1"] + member["p2"])))
            assert off <= POSITION_TOL + AT_TOLERANCE, "truss %s: witness point is %.3f mm from the web midpoint" % (
                truss.guid, off)
        else:
            off = distance_to_segment(point, member["p1"], member["p2"])
            assert off <= POSITION_TOL + AT_TOLERANCE, "truss %s: witness point is %.3f mm from the chord axis" % (
                truss.guid, off)


def named_demand_ratios(truss, member, witness, station, unit_station, load, critical_section):
    """Demand over limit of the check a witness names, for the carrier at `station` with gross load `load` (N).

    A chord point is read with its 1 mm tolerance, on both sides of a joint inside that range."""
    mode = witness["mode"]
    state = truss.state(station, witness["heavy_at"], unit_station)
    if member["kind"] == "web":
        assert mode in ("STRESS", "WEB_BUCKLING"), "truss %s: mode %s cannot name a web" % (truss.guid, mode)
        index = [w["guid"] for w in truss.webs].index(member["guid"])
        force = state.web_force(index, load)
        if mode == "STRESS":
            return [abs(force) / member["A"] / reference.STRESS_LIMIT]
        return [max(-force, 0.0) / float(truss.web_crit[index])]
    chord_index = 0 if member is truss.low else 1
    chord = truss.chords[chord_index]
    first, last = float(chord["x0"][0]), float(chord["x0"][-1] + chord["lengths"][-1])
    x = truss.station(np.array(witness["xyz_mm"], dtype=float))
    window = np.linspace(x - POSITION_TOL, x + POSITION_TOL, 41)
    joints = np.concatenate([chord["x0"], [last]])
    cuts = np.r_[joints, station, station + reference.WHEELBASE, unit_station, critical_section]
    window = np.concatenate([window, cuts[np.abs(cuts - x) <= POSITION_TOL + AT_TOLERANCE]])
    window = np.unique(np.clip(window, first, last))
    if mode == "STRESS":
        return [state.chord_stress(chord_index, float(v), load, side) / reference.STRESS_LIMIT
                for v in window for side in (1, -1)]
    assert mode == "SAG" and chord_index == 0, "truss %s: mode %s cannot name this chord" % (truss.guid, mode)
    return [state.sag(float(v), load) / truss.sag_limit for v in window]


def placement_ranges(truss, witness):
    """Feasible projections on the carrier station, with both reported stations allowed 1 mm."""
    h = witness.get("unit_s_mm")
    assert is_number(h), "truss %s: witness needs a numeric unit_s_mm" % truss.guid
    s = witness["s_mm"]
    tol = POSITION_TOL + AT_TOLERANCE
    a, b = max(0., h - tol), min(truss.span, h + tol)
    d, w = reference.CLEARANCE, reference.WHEELBASE
    intervals = []
    if a <= b:
        for lo, hi in ((a + d, truss.s_max), (a - w + d, b - d), (0., b - w - d)):
            lo, hi = max(0., s - tol, lo), min(truss.s_max, s + tol, hi)
            if lo <= hi:
                intervals.append((lo, hi))
    return intervals, a, b


@pytest.mark.parametrize("source", SOURCES)
def test_suspended_unit_clears_both_wheels(ws, source):
    """The witness places the suspended unit on the chord at least 400 mm from each wheel."""
    for truss, record in ws.pairs(source):
        witness = witness_of(truss, record)
        intervals, _, _ = placement_ranges(truss, witness)
        assert intervals, ("truss %s: unit at %.3f mm cannot clear both wheels at carrier station %.3f mm"
                           % (truss.guid, witness["unit_s_mm"], witness["s_mm"]))


def witness_placements(truss, witness, knots):
    """Test the reported configuration first, then the feasible edges of its station tolerances."""
    intervals, ha, hb = placement_ranges(truss, witness)
    s, h = witness["s_mm"], witness["unit_s_mm"]
    d, w = reference.CLEARANCE, reference.WHEELBASE
    visited = set()
    # An attained adverse configuration supplies additional interior probes when
    # they fall inside the reported tolerances. Other admissible witnesses have
    # exactly the same demand test.
    critical = {}
    for lo, hi in intervals:
        cuts = np.r_[np.clip(s, lo, hi), lo, hi, .5*(lo+hi), knots[(knots >= lo) & (knots <= hi)]]
        for z in cuts:
            if z not in critical:
                envelope = truss.envelope()
                upper, detail = envelope.evaluate(np.array([z]))
                code, section, unit, _ = detail[0]
                member, mode = envelope.labels[int(code)]
                matches = member["guid"] == witness["member"] and mode == witness["mode"]
                critical[z] = (upper[0], unit, section if matches else np.nan)
            capacity, adverse_unit, section = critical[z]
            for a, b in ((ha, min(hb, z-d)), (max(ha, z+d), min(hb, z+w-d)), (max(ha, z+w+d), hb)):
                if a <= b:
                    for unit in (np.clip(h, a, b), np.clip(adverse_unit, a, b), a, b, .5*(a+b)):
                        if (z, unit) not in visited:
                            visited.add((z, unit))
                            yield z, unit, capacity, section


@pytest.mark.parametrize("source", SOURCES)
def test_witness_demand_at_limit(ws, source):
    """Steel, unit and charted carrier load bring the named demand within 1.5% of its limit, or permanent load exceeds it."""
    for truss, record in ws.pairs(source):
        witness = witness_of(truss, record)
        member = witness_member(truss, witness)
        knots, ordinates = curve_arrays(truss, record)
        gap = math.inf
        for station, unit, capacity, section in witness_placements(truss, witness, knots):
            zero = capacity == 0.
            load = 0. if zero else 1000. * np.interp(station, knots, ordinates)
            ratios = np.asarray(named_demand_ratios(truss, member, witness, station, unit, load, section))
            miss = 1. - ratios.max() if zero else float(np.min(np.abs(ratios - 1.))) - WITNESS_TOL
            gap = min(gap, miss)
            if gap <= ROUND_OFF:
                break
        assert gap <= ROUND_OFF, (
            "truss %s: named %s demand at s=%.3f and unit_s_mm=%.3f misses the witness tolerance (gap %.5g)"
            % (truss.guid, witness["mode"], witness["s_mm"], witness["unit_s_mm"], gap))
