"""Compare the installed material and component axes independently of feature names."""

from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path

from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from solids import FreeCAD, Part
from frames import mounted_index


@lru_cache(maxsize=4)
def _solids(directory):
    shape = Part.read(str(Path(directory) / 'assembly.brep'))
    return shape.Solids


def _signature(solid):
    box = solid.optimalBoundingBox(False)
    return (box.XMin, box.YMin, box.ZMin, box.XMax, box.YMax, box.ZMax,
            *solid.CenterOfMass)


@lru_cache(maxsize=None)
def _pair_indices(expected_dir, actual_dir):
    left, right = _solids(expected_dir), _solids(actual_dir)
    if len(left) != len(right) or not left:
        return ()
    costs = cdist([_signature(s) for s in left], [_signature(s) for s in right],
                  metric='chebyshev')
    _, indices = linear_sum_assignment(costs)
    return tuple(map(int, indices))


def paired_solids(expected_dir, actual_dir):
    """Return comparator and mounted solids in the same geometric order."""
    order = _pair_indices(expected_dir, actual_dir)
    return _solids(expected_dir), [_solids(actual_dir)[i] for i in order]


def _datums(directory):
    return json.loads((Path(directory) / 'datums.json').read_text())


def _mounted_datums(expected_dir, actual_dir, order):
    def rows(folder):
        return json.loads((Path(folder) / "frames.json").read_text())
    expected, actual = rows(expected_dir), rows(actual_dir)
    worst, missing = 0., []
    for source in expected:
        if not source["leaf"] or not source["members"]:
            continue
        before = mounted_index(expected, source["members"])
        after = mounted_index(actual, [order[i] for i in source["members"]])
        if after is None:
            # A single, disconnected finished part still has one installation.
            # A union of independently mounted pieces does not supply its datum.
            missing.append(source["path"])
            continue
        def rotation(row):
            return FreeCAD.Rotation(*row["datum"][3:])
        angle = math.degrees((rotation(expected[before]).inverted() * rotation(actual[after])).Angle)
        worst = max(worst, min(angle, abs(360. - angle)))
    return worst, missing



def compare_inventory(expected_dir, actual_dir):
    """Measure occurrence inventory, validity and mounting-datum correspondence."""
    expected = _solids(expected_dir)
    actual = _solids(actual_dir)
    result = dict(count=len(actual), expected_count=len(expected), invalid=[],
                  angle_deg=None, missing_occurrences=[])
    result['invalid'] = [i for i, s in enumerate(actual)
                         if not s.isClosed() or not s.isValid() or s.Volume <= 0]
    result['invalid'] += _datums(actual_dir)['invalid']
    order = _pair_indices(expected_dir, actual_dir)
    if order:
        result['angle_deg'], result['missing_occurrences'] = _mounted_datums(expected_dir, actual_dir, order)
        if result['missing_occurrences']:
            result['angle_deg'] = math.inf
    return result


_certified_pairs = {}


def _local_id(solid):
    local = solid.copy()
    local.Placement = FreeCAD.Placement()
    identity = hashlib.sha256(local.exportBrepToString().encode()).digest()
    return local, identity


def _certified_deviation(a, b):
    """Reuse a boundary proof when only the relative mounting placement changes."""
    from boundary import TOL, boundary_error
    _, left = _local_id(a)
    local, right = _local_id(b)
    pose = a.Placement.inverse() * b.Placement
    previous = _certified_pairs.get((left, right))
    if previous is not None:
        old, error = previous
        box = local.BoundBox
        # Rigid displacement is convex on the box, so its largest possible
        # value is bounded by the largest displacement of the eight corners.
        vertices = product((box.XMin, box.XMax), (box.YMin, box.YMax), (box.ZMin, box.ZMax))
        drift = max((pose.multVec(FreeCAD.Vector(*p)) - old.multVec(FreeCAD.Vector(*p))).Length
                    for p in vertices)
        if error + drift <= TOL:
            return error + drift
    error = boundary_error(a, b)
    if error <= TOL:
        _certified_pairs[left, right] = pose, error
    return error


def _part_description(directory, index):
    rows = json.loads((Path(directory) / "frames.json").read_text())
    row = next(r for r in rows if r["leaf"] and index in r["members"])
    return f"{row['label']} ({'.'.join(row['path'])}, solid {row['members'].index(index) + 1})"


def _pair_description(expected_dir, actual_dir, i, j):
    ref = _part_description(expected_dir, i)
    got = _part_description(actual_dir, j)
    return f"comparator {ref} and release {got}"


def _witness_lower_bound(source, target):
    """Bound witness-to-skin distance below by distance to enclosing face boxes."""
    boxes = [face.BoundBox for face in target.Faces]
    points = [v.Point for v in source.Vertexes]
    points.extend(e.valueAt((e.FirstParameter + e.LastParameter) / 2) for e in source.Edges)
    largest = 0.
    for p in points:
        gap = min(math.hypot(max(box.XMin - p.x, p.x - box.XMax, 0.),
                             max(box.YMin - p.y, p.y - box.YMax, 0.),
                             max(box.ZMin - p.z, p.z - box.ZMax, 0.)) for box in boxes)
        largest = max(largest, gap)
    return largest


@lru_cache(maxsize=None)
def material_deviation(expected_dir, actual_dir):
    """Certify both skins and identify the installed pair on a failed comparison."""
    from boundary import TOL
    expected, actual = paired_solids(expected_dir, actual_dir)
    if not expected or len(actual) != len(expected):
        return math.inf, "No one-to-one installed solid correspondence"
    order = _pair_indices(expected_dir, actual_dir)
    upper = 0.
    for i, (a, b) in enumerate(zip(expected, actual)):
        try:
            bound = _certified_deviation(a, b)
        except ArithmeticError as exc:
            pair = _pair_description(expected_dir, actual_dir, i, order[i])
            raise ArithmeticError(f"{pair}: {exc}") from exc
        if bound > TOL:
            pair = _pair_description(expected_dir, actual_dir, i, order[i])
            # A disjoint tolerance ball sometimes proves only a strict threshold
            # crossing. Give a stronger physical bound when a real edge or
            # vertex lies outside the boxes enclosing all opposing faces.
            lower = max(_witness_lower_bound(a, b), _witness_lower_bound(b, a))
            message = f"{pair}: certified boundary deviation > {TOL} mm"
            if lower > TOL:
                message += f" (witness lower bound {lower:.9g} mm)"
            return max(bound, lower), message
        upper = max(upper, bound)
    return upper, ""
