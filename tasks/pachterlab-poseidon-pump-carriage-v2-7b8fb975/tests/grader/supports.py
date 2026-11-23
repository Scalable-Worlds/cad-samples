"""Measure the cylindrical supports, occupied stations and azimuths of the input."""

import io
import math
import os
from dataclasses import dataclass
from itertools import product

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import FreeCAD
import Part
import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepTools import BRepTools
from OCP.ShapeAnalysis import ShapeAnalysis_CanonicalRecognition
from OCP.TopoDS import TopoDS_Shape
from OCP.gp import gp_Cylinder


LENGTH_EQUALITY = 1e-5
ANGLE_EQUALITY = math.radians(1e-5)
LIMIT_EPS = 1e-6
TAU = math.tau


@dataclass
class SupportGroup:
    origin: np.ndarray
    axis: np.ndarray
    radius: float
    faces: list
    span: float = 0.0
    coverage: float = 0.0

    @property
    def eligible(self):
        return (1.25 - LIMIT_EPS <= self.radius <= 10 + LIMIT_EPS
                and self.span >= 1 - LIMIT_EPS and self.coverage >= math.pi - LIMIT_EPS)


def load_original(path):
    solid = Part.read(str(path))
    assert len(solid.Solids) == 1, f"Expected a single closed IGES solid in {path}"
    return solid.Solids[0]


def _inward_cylinder(face):
    # Canonical recognition checks the underlying surface, including NURBS.
    oc_face = TopoDS_Shape()
    BRepTools.Read_s(oc_face, io.BytesIO(face.exportBrepToString().encode()), BRep_Builder())
    cylinder = gp_Cylinder()
    if not ShapeAnalysis_CanonicalRecognition(oc_face).IsCylinder(LENGTH_EQUALITY, cylinder):
        return None
    origin = np.array(cylinder.Location().Coord())
    axis = np.array(cylinder.Axis().Direction().Coord())
    radius = cylinder.Radius()
    u0, u1, v0, v1 = face.ParameterRange
    uv = ((u0 + u1) / 2, (v0 + v1) / 2)
    from_axis = np.asarray(face.valueAt(*uv)) - origin
    radial = from_axis - axis * (from_axis @ axis)
    # normalAt includes the orientation of the face in the original solid.
    if radial @ np.asarray(face.normalAt(*uv)) >= 0:
        return None
    return origin, axis, radius


def _coincident(group, origin, axis, radius, box_corners):
    if (abs(group.radius - radius) > LENGTH_EQUALITY
            or np.linalg.norm(np.cross(group.axis, axis)) > math.sin(ANGLE_EQUALITY)):
        return False
    stations = (box_corners - group.origin) @ group.axis
    ends = group.origin + np.outer([stations.min(), stations.max()], group.axis)
    return np.linalg.norm(np.cross(ends - origin, axis), axis=1).max() <= LENGTH_EQUALITY


def _measure_footprint(group, faces):
    # The in-plane basis may have any phase: interval union on a circle is invariant.
    axis = group.axis
    first = np.eye(3)[np.argmin(abs(axis))]
    first = first - (first @ axis) * axis
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    intervals = []
    low, high = math.inf, -math.inf
    for index in group.faces:
        vertices, triangles = faces[index].tessellate(0.01)
        xyz = np.asarray([tuple(p) for p in vertices]) - group.origin
        tri = np.asarray(triangles)
        stations = xyz[np.unique(tri)] @ axis
        low, high = min(low, stations.min()), max(high, stations.max())
        angles = np.mod(np.arctan2(xyz @ second, xyz @ first), TAU)
        points = np.sort(angles[tri], axis=1)
        # Each triangle occupies the complement of its widest circular gap.
        gaps = np.roll(points, -1, axis=1) - points
        gaps[:, 2] += TAU
        widest = np.argmax(gaps, axis=1)
        rows = np.arange(len(points))
        begins = points[rows, (widest + 1) % 3]
        ends = points[rows, widest]
        wraps = ends < begins
        intervals.append(np.column_stack((begins, np.where(wraps, TAU, ends))))
        intervals.append(np.column_stack((np.zeros(wraps.sum()), ends[wraps])))
    segments = np.concatenate(intervals)
    segments = segments[np.argsort(segments[:, 0])]
    # Only the portion past all preceding right endpoints adds new coverage.
    previous_end = np.r_[0.0, np.maximum.accumulate(segments[:-1, 1])]
    total = np.maximum(0.0, segments[:, 1] - np.maximum(previous_end, segments[:, 0])).sum()
    group.span, group.coverage = high - low, total


def groups_for(shape):
    """Return every inward cylindrical support, with zero-based face indices."""
    box = shape.BoundBox
    corners = np.array(list(product((box.XMin, box.XMax), (box.YMin, box.YMax),
                                    (box.ZMin, box.ZMax))))
    faces = shape.Faces
    groups = []
    for index, face in enumerate(faces):
        cylinder = _inward_cylinder(face)
        if cylinder is None:
            continue
        origin, axis, radius = cylinder
        for group in groups:
            if _coincident(group, origin, axis, radius, corners):
                group.faces.append(index)
                break
        else:
            groups.append(SupportGroup(origin, axis, radius, [index]))
    for group in groups:
        _measure_footprint(group, faces)
    return groups
