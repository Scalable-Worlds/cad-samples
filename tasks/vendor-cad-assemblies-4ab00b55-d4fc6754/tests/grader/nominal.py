"""Develop the selected skin of a registered analytic sheet-metal solid."""
import os
from collections import defaultdict, deque
from functools import lru_cache

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import FreeCAD
import Part
import numpy as np
from scipy.spatial import ConvexHull


def coordinates(vector):
    return np.array(tuple(vector))


def outward(face, point):
    parameters = face.Surface.parameter(FreeCAD.Vector(*point))
    return coordinates(face.normalAt(*parameters))


def flat_map(origin, axis, transverse, start, along, across):
    matrix = np.outer(along, axis) + np.outer(across, transverse)

    def map_points(points):
        return (points - origin) @ matrix.T + start

    return map_points


def bend_map(center, axis, radial_start, radial_finish, radius, start, along, across):
    def map_points(points):
        offsets = points - center
        angles = np.arctan2(offsets @ radial_finish, offsets @ radial_start)
        return start + np.outer(offsets @ axis, along) + radius * np.outer(angles, across)

    return map_points


def radial_direction(point, center, axis):
    offset = point - center
    radial = offset - axis * (offset @ axis)
    return radial / np.linalg.norm(radial)


@lru_cache(maxsize=None)
def registered_stock(master_path, stl_path):
    record = np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")])
    vertices = np.fromfile(stl_path, record, offset=84)["vertices"]
    coordinate_limit = 500.0
    largest_coordinate = float(np.max(np.abs(vertices)))
    print(f"METRIC name=stock_abs_coordinate_mm value={largest_coordinate} bar={coordinate_limit} sense=max",
          flush=True)
    assert largest_coordinate <= coordinate_limit, f"Stock placement exceeds ±{coordinate_limit:g} mm"
    mesh_min = vertices.min(axis=(0, 1))
    solid = Part.read(str(master_path))
    box = solid.BoundBox
    translation = mesh_min - np.array([box.XMin, box.YMin, box.ZMin])
    faces = solid.Faces
    incident = defaultdict(list)
    for index, face in enumerate(faces):
        for edge in face.Edges:
            incident[solid.findSubShape(edge)[1]].append((index, edge))

    tangent = defaultdict(dict)
    for neighbors in incident.values():
        if len(neighbors) != 2:
            continue
        (i, edge), (j, _) = neighbors
        surfaces = (faces[i].Surface, faces[j].Surface)
        if not any(isinstance(s, Part.Cylinder) for s in surfaces):
            continue
        if not any(isinstance(s, Part.Plane) for s in surfaces):
            continue
        point = coordinates(edge.valueAt((edge.FirstParameter + edge.LastParameter) / 2))
        if outward(faces[i], point) @ outward(faces[j], point) > 1 - 1e-9:
            tangent[i][j] = point
            tangent[j][i] = point
    return solid, translation, faces, incident, tangent


def skin_boundaries(solid, incident, maps):
    boundaries = []
    at_vertex = defaultdict(list)
    for neighbors in incident.values():
        skin_faces = [(i, e) for i, e in neighbors if i in maps]
        if len(skin_faces) != 1:
            continue
        index, edge = skin_faces[0]
        if isinstance(edge.Curve, Part.Line):
            sampled = [edge.valueAt(edge.FirstParameter), edge.valueAt(edge.LastParameter)]
        else:
            sampled = edge.discretize(Deflection=0.0004)
        points = np.array([coordinates(p) for p in sampled])
        # Edges in opposite face orientations can expose vertices in either order.
        # Chain by CAD vertex identity so seams introduce no polygon snapping.
        vertices = edge.Vertexes
        ends = []
        for point in (points[0], points[-1]):
            vertex = min(vertices, key=lambda v: np.linalg.norm(coordinates(v.Point) - point))
            ends.append(solid.findSubShape(vertex)[1])
        for vertex in ends:
            at_vertex[vertex].append(len(boundaries))
        boundaries.append((*ends, maps[index](points)))

    rings = []
    available = set(range(len(boundaries)))
    while available:
        index = available.pop()
        start, vertex, xy = boundaries[index]
        portions = [xy[:-1]]
        while vertex != start:
            index, = (k for k in at_vertex[vertex] if k in available)
            available.remove(index)
            first, last, xy = boundaries[index]
            if last == vertex:
                xy = xy[::-1]
                first, last = last, first
            portions.append(xy[:-1])
            vertex = last
        rings.append(np.concatenate(portions).tolist())
    return rings


def develop(master_path, stl_path, datum):
    """Return neutral-space trim loops and signed mid-angle generators in the drawing frame."""
    solid, translation, faces, incident, tangent = registered_stock(master_path, stl_path)
    pick = np.asarray(datum) - translation
    vertex = Part.Vertex(FreeCAD.Vector(*pick))
    planes = [i for i, face in enumerate(faces) if isinstance(face.Surface, Part.Plane)]
    root = min(planes, key=lambda i: faces[i].distToShape(vertex)[0])
    face = faces[root]
    normal = outward(face, coordinates(face.CenterOfMass))
    origin = pick - normal * ((pick - coordinates(face.CenterOfMass)) @ normal)
    first = np.eye(3)[np.argmin(np.round(np.abs(normal), 6))]
    first = first - (first @ normal) * normal
    first /= np.linalg.norm(first)
    maps = {root: flat_map(origin, first, np.cross(normal, first), np.zeros(2), [1, 0], [0, 1])}
    bend_count = {root: 0}
    pending = deque([root])
    brakes = {"INNER": [], "OUTER": []}
    while pending:
        current = pending.popleft()
        for bend, entry in tangent[current].items():
            if bend in maps:
                continue
            other, = (i for i in tangent[bend] if i != current)
            patch = faces[bend]
            axis = coordinates(patch.Surface.Axis)
            center = coordinates(patch.Surface.Center)
            r = patch.Surface.Radius
            d0 = radial_direction(entry, center, axis)
            d1 = radial_direction(tangent[bend][other], center, axis)
            sign = np.sign(d0 @ outward(patch, entry))
            neutral = r - sign * 0.75

            # For a quarter-turn the entry tangent is the exit radial direction.
            # The two lines at the tangencies share the same axial origin.
            entry_origin = center + r * d0
            xy = maps[current](entry_origin)
            along = maps[current](entry_origin + axis) - xy
            across = maps[current](entry_origin + d1) - xy
            maps[bend] = bend_map(center, axis, d0, d1, neutral, xy, along, across)
            if other not in maps:
                finish = xy + across * (np.pi / 2) * neutral
                maps[other] = flat_map(center + r * d1, axis, -d0, finish, along, across)
                bend_count[other] = bend_count[current] + 1
                pending.append(other)

            mid = center + r * (d0 + d1) / np.sqrt(2)
            extents = [(coordinates(v.Point) - mid) @ axis for v in patch.Vertexes]
            generator = Part.makeLine(FreeCAD.Vector(*(mid + (min(extents) - 1) * axis)),
                                      FreeCAD.Vector(*(mid + (max(extents) + 1) * axis)))
            for edge in patch.common(generator).Edges:
                endpoints = np.array([coordinates(v.Point) for v in edge.Vertexes])
                brakes["OUTER" if sign > 0 else "INNER"].append(maps[bend](endpoints).tolist())
    # Check the scope of the job on the nominal stock, independently of the DXF.
    # Breadth-first traversal gives the shortest bend path to each panel.
    depth = max(bend_count.values())
    print(f"METRIC name=stock_bend_depth value={depth} bar=3 sense=max", flush=True)
    assert depth <= 3, "Nominal panel requires more than three bends from the datum"
    rings = skin_boundaries(solid, incident, maps)
    points = np.concatenate(rings)
    hull = points[ConvexHull(points).vertices]
    diameter = max(np.linalg.norm(hull - p, axis=1).max() for p in hull)
    print(f"METRIC name=stock_developed_diameter_mm value={diameter} bar=700 sense=max", flush=True)
    assert diameter <= 700, "Nominal developed blank exceeds 700 mm between points"
    return {"CUT": rings, **brakes}
