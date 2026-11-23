"""Use the design cylinders to identify and section their corresponding STL walls."""

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import FreeCAD
import Part
import numpy as np


@dataclass
class Land:
    c: np.ndarray
    d: np.ndarray
    r: float
    triangles: np.ndarray


@lru_cache(maxsize=1)
def _design_cylinders():
    """Read the analytic cylinder faces extracted from each matching design solid."""
    result = {}
    for path in sorted((Path(__file__).parent / "cylinders").glob("*.brep")):
        groups = {}
        for face in Part.read(str(path)).Faces:
            surface = face.Surface
            radius = surface.Radius
            if not 3.5 <= 2 * radius <= 40:
                continue

            axis = np.array(surface.Axis)
            axis /= np.linalg.norm(axis)
            center = np.array(surface.Center)
            center -= axis * (axis @ center)
            u0, u1, v0, v1 = face.ParameterRange
            u, v = (u0 + u1) / 2, (v0 + v1) / 2
            radial = np.array(face.valueAt(u, v)) - center
            if np.array(face.normalAt(u, v)) @ radial >= 0:
                continue

            if axis[np.argmax(np.abs(axis))] < 0:
                axis = -axis
            # Repeated STEP surfaces can use opposite parameterizations and
            # different origins along the same straight line.
            key = tuple(np.r_[center, axis, radius].round(7))
            if key not in groups:
                groups[key] = (center, axis, radius, [])
            groups[key][3].append(face)
        result[path.stem] = list(groups.values())
    return result


def _merge(intervals, roundoff=0):
    merged = []
    for lo, hi in sorted(intervals):
        if hi <= lo:
            continue
        if merged and lo <= merged[-1][1] + roundoff:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([float(lo), float(hi)])
    return merged


def _angular_measure(radial_vectors, axis):
    # Every facet spans less than half a turn. Split at the zero meridian
    # before taking the union so disjoint axial portions contribute only once.
    e = np.eye(3)[np.argmin(np.abs(axis))]
    e -= (e @ axis) * axis
    e /= np.linalg.norm(e)
    angles = np.unwrap(
        np.arctan2(radial_vectors @ np.cross(axis, e), radial_vectors @ e),
        axis=1,
    )
    spans = []
    turn = 2 * np.pi
    for low, high in zip(angles.min(axis=1), angles.max(axis=1)):
        start = low % turn
        end = start + (high - low)
        spans.append((start, min(turn, end)))
        if end > turn:
            spans.append((0, end - turn))
    return sum(high - low for low, high in _merge(spans))


def get_lands(key, original_triangles, R, t, placed_triangles=None):
    """Find qualifying walls on a print and place them in its survey frame."""
    tri = np.asarray(original_triangles, dtype=float)
    placed = tri @ R.T + t if placed_triangles is None else np.asarray(placed_triangles, dtype=float)
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    points = tri.mean(axis=1)
    lands = []
    for center, axis, radius, trims in _design_cylinders()[key]:
        offset = tri - center
        radial = offset - (offset @ axis)[:, :, None] * axis
        distances = np.linalg.norm(radial, axis=2)
        mid_radial = radial.mean(axis=1)
        inward = (normals * mid_radial).sum(axis=1) < 0
        # Original STL vertices are float32, including those at CAD trim
        # intersections. The chordal envelope bounds the distance to the
        # design cylinder. The exact trim determines which face owns the facet.
        precision = np.finfo(np.float32).eps * max(1.0, np.abs(tri).max())
        in_envelope = (np.abs(distances - radius) <= 0.05 + precision).all(axis=1)
        # A shoulder triangle can share all three of its vertices with a
        # cylindrical rim. Its chords can span the interior of the bore.
        neighbor = np.roll(radial, -1, axis=1)
        radial_dot = (radial * neighbor).sum(axis=2)
        radii_product = distances * np.roll(distances, -1, axis=1)
        within_span = (radial_dot >= np.cos(np.deg2rad(16)) * radii_product).all(axis=1)
        face_ids = []
        for i in np.flatnonzero(inward & in_envelope & within_span):
            radial_direction = mid_radial[i] / np.linalg.norm(mid_radial[i])
            on_surface = center + (axis @ points[i]) * axis + radius * radial_direction
            location = FreeCAD.Vector(*on_surface)
            if any(face.isPartOfDomain(*face.Surface.parameter(location)) for face in trims):
                face_ids.append(i)
        if not face_ids or np.degrees(_angular_measure(radial[face_ids], axis)) < 160:
            continue
        d = R @ axis
        p = R @ center + t
        c = p - (p @ d) * d
        lands.append(Land(c, d, radius, placed[face_ids]))
    return lands


def wall_rows(land, d, normals):
    """Intersect the actual land triangles with each outward radial half-plane."""
    offset = land.triangles - land.c
    heights = land.triangles @ d
    next_height = np.roll(heights, -1, axis=1)
    mid_radial = offset.mean(axis=1)
    arithmetic_error = 32 * np.finfo(float).eps * max(1.0, np.max(np.abs(land.triangles)))
    rows = []
    for direction in normals:
        signed = offset @ np.cross(d, direction)
        signed[np.abs(signed) <= arithmetic_error] = 0
        adjacent = np.roll(signed, -1, axis=1)
        crossed = ((signed < 0) & (adjacent > 0)) | ((signed > 0) & (adjacent < 0))
        fraction = np.divide(signed, signed - adjacent, where=crossed, out=np.zeros_like(signed))
        cut = heights + fraction * (next_height - heights)

        # Include any coincident vertex or axial edge. Reject intersections
        # with the other half of the cutting plane.
        low = np.minimum(
            np.min(np.where(crossed, cut, np.inf), axis=1),
            np.min(np.where(signed == 0, heights, np.inf), axis=1),
        )
        high = np.maximum(
            np.max(np.where(crossed, cut, -np.inf), axis=1),
            np.max(np.where(signed == 0, heights, -np.inf), axis=1),
        )
        hit = (mid_radial @ direction > 0) & (high > low)
        rows.append(_merge(zip(low[hit], high[hit]), arithmetic_error))
    return rows
