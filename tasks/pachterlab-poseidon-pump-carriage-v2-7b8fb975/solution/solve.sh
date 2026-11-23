#!/bin/bash
set -euo pipefail

mkdir -p /app/out
rm -f -- /app/out/relieve.py /app/out/fit.step
cat > /app/out/relieve.py <<'PYTHON_RELIEF_SOURCE'
#!/opt/conda/bin/python
"""Apply a 0.20 mm radial print allowance to qualifying trimmed IGES walls."""
import argparse
from dataclasses import dataclass, field
import io
import math
import os
from pathlib import Path
import tempfile

import FreeCAD
import Part
import numpy as np
from OCP.BRep import BRep_Builder
from OCP.BRepClass import BRepClass_FaceClassifier
from OCP.BRepTools import BRepTools
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.IntCurveSurface import IntCurveSurface_Tangent
from OCP.TopAbs import TopAbs_IN
from OCP.TopoDS import TopoDS, TopoDS_Shape
from OCP.gp import gp_Dir, gp_Lin, gp_Pnt, gp_Pnt2d
from scipy.optimize import minimize_scalar
from scipy.stats import qmc


# The final term only covers floating-point roundoff at a comparison.
_INCLUSION_EPS = 1e-6 + 1e-12


def _unit_np(vector):
    v = np.asarray(vector, dtype=float)
    return v / np.linalg.norm(v)


def _canonical_axis(axis):
    axis = _unit_np(axis)
    return axis if axis[np.argmax(np.abs(axis))] > 0 else -axis


def _param_sites(lo, hi, knots, degree):
    # Stay strictly inside each polynomial span to avoid one-sided
    # derivatives or singularities at the boundaries of a carrier.
    breaks = [lo] + [k for k in knots if lo < k < hi] + [hi]
    fractions = (np.arange(max(3, 2 * degree + 1)) + 0.5) / max(3, 2 * degree + 1)
    return np.concatenate([a + (b - a) * fractions
                           for a, b in zip(breaks, breaks[1:]) if b > a])


def _inward_cylinder(face):
    """Return a point on the axis, axis direction, and radius in world coordinates."""
    surf = face.Surface
    u0, u1, v0, v1 = face.ParameterRange
    if surf.TypeId == 'Part::GeomCylinder':
        axis = _canonical_axis(surf.Axis)
        center = np.asarray(surf.Center, dtype=float)
        u, v = (u0 + u1) / 2, (v0 + v1) / 2
        radial = np.asarray(face.valueAt(u, v), dtype=float) - center
        radial -= axis * np.dot(radial, axis)
        if np.dot(radial, np.asarray(face.normalAt(u, v))) < 0:
            return center, axis, surf.Radius
        return None
    if surf.TypeId in ('Part::GeomPlane', 'Part::GeomSphere', 'Part::GeomToroid'):
        return None
    if surf.isPlanar(1e-9):
        return None

    # For a right circular cylinder all oriented unit normals lie in one
    # vector plane. After projecting normal to its axis, p = c +/- r*n.
    # This also works when neither surface parameter is an axial coordinate.
    if hasattr(surf, 'getUKnots'):
        us = _param_sites(u0, u1, surf.getUKnots(), surf.UDegree)
        vs = _param_sites(v0, v1, surf.getVKnots(), surf.VDegree)
    else:
        us = _param_sites(u0, u1, [], 4)
        vs = _param_sites(v0, v1, [], 4)
    sites = [(u, v) for u in us for v in vs if face.isPartOfDomain(u, v)]
    if len(sites) < 9:
        # Very narrow or diagonal trims can miss most of a rectangular UV grid.
        # Interpolate within triangles already meshed on the occupied domain.
        _, triangles = face.tessellate(0.01)
        if triangles:
            nodes = np.asarray(face.getUVNodes())
            step = max(1, len(triangles) // 96)
            extra = nodes[np.asarray(triangles[::step])].mean(axis=1)
            sites.extend((u, v) for u, v in extra if face.isPartOfDomain(u, v))
    samples, normals = [], []
    for u, v in sites:
        try:
            normal = np.asarray(face.normalAt(u, v), dtype=float)
            point = np.asarray(face.valueAt(u, v), dtype=float)
        except (RuntimeError, ValueError):
            continue
        if np.all(np.isfinite(point)) and np.all(np.isfinite(normal)) and np.linalg.norm(normal) > 0.99:
            samples.append(point)
            normals.append(normal)
    if len(normals) < 3:
        return None
    p, n = np.asarray(samples), np.asarray(normals)
    _, singular, vt = np.linalg.svd(n, full_matrices=False)
    if singular[1] < 1e-11 * singular[0]:
        return None
    axis = _canonical_axis(vt[-1])
    if np.max(np.abs(n @ axis)) > 1.75e-7:
        return None
    pmean, nmean = p.mean(axis=0), n.mean(axis=0)
    dp = p - pmean
    dp -= np.outer(dp @ axis, axis)
    dn = n - nmean
    radius = -np.sum(dp * dn) / np.sum(dn * dn)
    if not (1.25 - _INCLUSION_EPS <= radius <= 10 + _INCLUSION_EPS):
        return None
    if np.max(np.linalg.norm(dp + radius * dn, axis=1)) > 1e-5:
        return None
    return pmean + radius * nmean, axis, radius


@dataclass
class CylinderGroup:
    origin: np.ndarray
    axis: np.ndarray
    radius: float
    indices: list = field(default_factory=list)
    span: float = 0.0
    coverage: float = 0.0


def _same_support(group, center, axis, radius, solid, box_center, box_halfsize):
    if abs(radius - group.radius) > 1e-5:
        return False
    sine = np.linalg.norm(np.cross(group.axis, axis))
    if sine > math.sin(math.radians(1e-5)):
        return False
    # Compare the lines at both ends of the part, independently of which
    # point each carrier uses as its origin.
    midpoint = group.origin + group.axis * np.dot(box_center - group.origin, group.axis)
    extent = np.dot(np.abs(group.axis), box_halfsize)
    def transverse(point):
        return np.linalg.norm(np.cross(point - center, axis))

    if transverse(midpoint) > 1e-5 + extent * sine:
        return False
    if max(transverse(midpoint - extent * group.axis),
           transverse(midpoint + extent * group.axis)) <= 1e-5:
        return True
    # For nearly coincident, slightly tilted supports the corners of the
    # bounding box may extend beyond the material. Use actual axial extrema
    # to resolve that borderline case.
    rotation = FreeCAD.Rotation(FreeCAD.Vector(*group.axis), FreeCAD.Vector(0, 0, 1))
    axial_box = solid.transformed(rotation.toMatrix()).optimalBoundingBox(False, False)
    station = np.dot(midpoint, group.axis)
    return all(transverse(midpoint + (z - station) * group.axis) <= 1e-5
               for z in (axial_box.ZMin, axial_box.ZMax))


def _wall_extent(group, faces):
    axis = group.axis
    ref = np.eye(3)[np.argmin(np.abs(axis))]
    xdir = _unit_np(np.cross(axis, ref))
    ydir = np.cross(axis, xdir)
    tau = 2 * math.pi
    intervals = []
    zmin, zmax = math.inf, -math.inf
    for index in group.indices:
        vertices, facets = faces[index].tessellate(0.01)
        if not facets:
            continue
        p = np.array(vertices, dtype=float) - group.origin
        facets = np.asarray(facets, dtype=int)
        z = (p @ axis)[facets]
        zmin, zmax = min(zmin, z.min()), max(zmax, z.max())
        angles = np.mod(np.arctan2(p @ ydir, p @ xdir), tau)
        triples = np.sort(angles[facets], axis=1)
        # The complement of the largest gap is the short covering arc.
        end = np.column_stack((triples[:, 1:], triples[:, 0] + tau))
        largest = np.argmax(end - triples, axis=1)
        rows = np.arange(len(facets))
        start = triples[rows, (largest + 1) % 3]
        stop = triples[rows, largest]
        wrap = stop < start
        intervals.extend(zip(start[~wrap], stop[~wrap]))
        intervals.extend(zip(start[wrap], np.full(wrap.sum(), tau)))
        intervals.extend(zip(np.zeros(wrap.sum()), stop[wrap]))
    if not intervals:
        return
    length, right = 0.0, 0.0
    for lo, hi in sorted(intervals):
        if hi > right:
            length += hi - max(right, lo)
            right = hi
    group.coverage, group.span = length, zmax - zmin


def find_groups(solid):
    faces = solid.Faces
    box = solid.BoundBox
    box_center = np.asarray(box.Center, dtype=float)
    box_halfsize = np.array([box.XLength, box.YLength, box.ZLength]) / 2
    groups = []
    for index, face in enumerate(faces):
        support = _inward_cylinder(face)
        if support is None:
            continue
        center, axis, radius = support
        if not (1.25 - _INCLUSION_EPS <= radius <= 10 + _INCLUSION_EPS):
            continue
        group = next((g for g in groups if _same_support(g, center, axis, radius, solid,
                                                          box_center, box_halfsize)), None)
        if group is None:
            group = CylinderGroup(center + axis * np.dot(box_center - center, axis), axis, radius)
            groups.append(group)
        group.indices.append(index)
    qualified = []
    for group in groups:
        _wall_extent(group, faces)
        if group.span >= 1.0 - _INCLUSION_EPS and group.coverage >= math.pi - _INCLUSION_EPS:
            qualified.append(group)
    return qualified


def make_cutters(solid, group):
    """Sweep each trimmed patch in the direction of increasing cylindrical radius.

    At every station the target wall is a homothety in the plane normal to
    the axis. A single 3D affine map therefore acts on all poles, retaining
    their original rational weights and the complete trim topology.
    """
    axis = np.asarray(group.axis, dtype=float)
    origin = np.asarray(group.origin, dtype=float)
    scale_delta = (0.20 / group.radius) * (np.eye(3) - np.outer(axis, axis))
    affine = np.eye(4)
    affine[:3, :3] += scale_delta
    affine[:3, 3] = -scale_delta @ origin
    transform = FreeCAD.Matrix(*affine.ravel())
    cutters = []
    faces = solid.Faces
    for index in group.indices:
        # Give the near and far trim curves the same spline parameterization
        # before applying the affine map, also when the input is analytic.
        near = faces[index].copy().toNurbs()
        far = near.transformGeometry(transform)
        boundaries = [near, far]
        occurrences = [edge for wire in near.Wires for edge in wire.OrderedEdges]
        for edge, outer_edge in zip(near.Edges, far.Edges, strict=True):
            # The wire traverses an actual periodic seam in both directions.
            # A one-sided trim can also lie on the carrier seam and needs a side.
            if edge.Degenerated or sum(edge.isSame(e) for e in occurrences) != 1:
                continue
            # Matching points on these two curves have the same (azimuth, z).
            boundaries.append(Part.makeRuledSurface(edge, outer_edge))
        shell = Part.makeShell(boundaries)
        if not (shell.isClosed() and shell.isValid()):
            sewn = Part.makeCompound(boundaries)
            sewn.sewShape()
            if len(sewn.Shells) != 1:
                raise RuntimeError(f"Could not sew the radial sweep of face {index + 1}")
            shell = sewn.Shells[0]
        cutter = Part.makeSolid(shell)
        if cutter.Volume < 0:
            cutter.reverse()
        if not (cutter.isClosed() and cutter.isValid() and cutter.Volume > 0):
            raise RuntimeError(f"Invalid radial sweep on face {index + 1}")
        cutters.append(cutter)
    return cutters


# Independent material checks use the original trim geometry and the imported STEP.



def _occ_brep(shape):
    result = TopoDS_Shape()
    stream = io.BytesIO(shape.exportBrepToString().encode("ascii"))
    BRepTools.Read_s(result, stream, BRep_Builder())
    return result


class _ParitySolid:
    directions = tuple(gp_Dir(*d) for d in (
        (0.33585, 0.7985, 0.50049), (0.59211, -0.18258, 0.78342),
        (-0.56118, 0.61272, 0.56143),
    ))

    def __init__(self, solid):
        self.shape = _occ_brep(solid)
        box = solid.BoundBox
        self.center = np.array(box.Center)
        self.radius = box.DiagonalLength * 0.5 + 1.0
        self.rays = None
        if solid.countElement("Face") < 600:
            self.rays = IntCurvesFace_ShapeIntersector()
            self.rays.Load(self.shape, 1e-7)
            return

        # Spatially partition very large STEP solids. Every face in a leaf
        # intersected by the ray still goes through the same OCCT intersector.
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        from OCP.TopAbs import TopAbs_FACE
        from OCP.TopExp import TopExp_Explorer
        from OCP.TopoDS import TopoDS_Compound
        faces, bounds = [], []
        explorer = TopExp_Explorer(self.shape, TopAbs_FACE)
        while explorer.More():
            face = explorer.Current()
            b = Bnd_Box()
            BRepBndLib.Add_s(face, b, False)
            faces.append(face)
            bounds.append(b.Get())
            explorer.Next()
        bounds = np.asarray(bounds)
        centers = (bounds[:, :3] + bounds[:, 3:]) / 2
        chunks = []
        leaves = []
        builder = BRep_Builder()

        def partition(ids):
            lo = bounds[ids, :3].min(axis=0)
            hi = bounds[ids, 3:].max(axis=0)
            if len(ids) <= 16:
                compound = TopoDS_Compound()
                builder.MakeCompound(compound)
                for i in ids:
                    builder.Add(compound, faces[i])
                chunks.append(compound)
                leaves.append((lo - 1e-6, hi + 1e-6))
            else:
                dim = np.argmax(hi - lo)
                ids = ids[np.argsort(centers[ids, dim])]
                mid = len(ids) // 2
                partition(ids[:mid])
                partition(ids[mid:])
        partition(np.arange(len(faces)))
        self.boxes = np.asarray(leaves)
        self.chunks = chunks
        self.intersectors = {}

    def _candidates(self, xyz, direction, reach):
        if self.rays is not None:
            yield self.rays
            return
        parameters = (self.boxes - xyz) / np.array(direction.Coord())
        entry = parameters.min(axis=1).max(axis=1)
        leave = parameters.max(axis=1).min(axis=1)
        for i in np.flatnonzero(np.maximum(entry, 0) <= np.minimum(leave, reach)):
            if i not in self.intersectors:
                rays = IntCurvesFace_ShapeIntersector()
                rays.Load(self.chunks[i], 1e-7)
                self.intersectors[i] = rays
            yield self.intersectors[i]

    def contains(self, point):
        start = gp_Pnt(*point)
        reach = self.radius + np.linalg.norm(point - self.center)
        odd = 0
        for direction in self.directions:
            hits = []
            line = gp_Lin(start, direction)
            for rays in self._candidates(point, direction, reach):
                rays.Perform(line, 0.0, reach)
                if not rays.IsDone():
                    raise RuntimeError("OCCT ray intersection failed")
                hits.extend(rays.WParameter(i) for i in range(1, rays.NbPnt() + 1)
                            if rays.Transition(i) != IntCurveSurface_Tangent
                            and rays.WParameter(i) > 0.0)
            crossings, previous = 0, -math.inf
            for distance in sorted(hits):
                if distance - previous > 1e-6:
                    crossings += 1
                previous = distance
            odd += crossings % 2
        return odd >= 2


def _occ_trim_edges(face):
    # A trim can end at the surface's parametric seam and still bound a patch.
    uses = [edge for wire in face.Wires for edge in wire.OrderedEdges]
    return [edge for edge in face.Edges if not edge.Degenerated
            and not (edge.isSeam(face) and sum(e.isSame(edge) for e in uses) >= 2)]


def _occ_sites(face, count):
    occ_face = TopoDS.Face(_occ_brep(face))
    u0, u1, v0, v1 = face.ParameterRange
    samples = qmc.Halton(2, scramble=False).random(32 * count + 2)[1:]
    found = 0
    for a, b in samples:
        u, v = u0 + a * (u1 - u0), v0 + b * (v1 - v0)
        if BRepClass_FaceClassifier(occ_face, gp_Pnt2d(u, v), 1e-9).State() == TopAbs_IN:
            yield u, v
            found += 1
            if found == count:
                return


def _occ_frame(group, point):
    delta = point - group.origin
    base = group.origin + group.axis * np.dot(delta, group.axis)
    v = point - base
    radius = np.linalg.norm(v)
    return base, v / radius if radius else v, radius


def _occ_probes(faces, groups):
    group_for = {i: g for g in groups for i in g.indices}
    for i, face in enumerate(faces):
        g = group_for.get(i)
        for j, (u, v) in enumerate(_occ_sites(face, 8)):
            p = np.array(face.valueAt(u, v))
            n = np.array(face.normalAt(u, v))
            for step in (-0.07, 0.07):
                yield p + step*n, (i + 1, "normal", j + 1, step)
        if g is None:
            continue
        for j, (u, v) in enumerate(_occ_sites(face, 16)):
            c, d, _ = _occ_frame(g, np.array(face.valueAt(u, v)))
            for step in (-0.07, 0.06, 0.12, 0.16, 0.24, 0.28):
                yield c + (g.radius + step)*d, (i + 1, "radial", j + 1, step)
        for j, edge in enumerate(_occ_trim_edges(face)):
            low, high = edge.ParameterRange
            for fraction in (0.15, 0.50, 0.85):
                p = np.array(edge.valueAt(low + fraction*(high - low)))
                c, d, _ = _occ_frame(g, p)
                tangent = np.cross(g.axis, d)
                for step in (-0.12, 0.12):
                    angle = step / g.radius
                    rotated = math.cos(angle)*d + math.sin(angle)*tangent
                    for dr in (0.10, 0.25):
                        r = g.radius + dr
                        key = (i + 1, j + 1, fraction, step, dr)
                        yield c + step*g.axis + r*d, ("axial edge", *key)
                        yield c + r*rotated, ("hoop edge", *key)


def _occ_contains_projection(faces, g, projected):
    vertex = Part.Vertex(FreeCAD.Vector(*projected))
    return any(vertex.distToShape(faces[i])[0] <= 1e-5 for i in g.indices)


def _occ_relief_contains(faces, groups, point):
    for g in groups:
        c, d, r = _occ_frame(g, point)
        if g.radius < r < g.radius + 0.20:
            if _occ_contains_projection(faces, g, c + g.radius*d):
                return True
    return False


def _occ_edge_boundary(faces, group):
    edges = []
    for i in group.indices:
        for edge in _occ_trim_edges(faces[i]):
            for j, (other, n) in enumerate(edges):
                if edge.isSame(other):
                    edges[j] = (other, n + 1)
                    break
            else:
                edges.append((edge, 1))
    return [e for e, uses in edges if uses == 1]


def _occ_near_side(point, group, edge, tolerance):
    # At each curve parameter the side is a radial line segment. The closest
    # point on that segment follows from its scalar projection. Minimize the
    # residual over the trim curve, including both ends and any spline knots.
    base, direction, _ = _occ_frame(group, point)
    projected = base + group.radius*direction
    dist, pairs, _ = edge.distToShape(Part.Vertex(FreeCAD.Vector(*projected)))
    if dist > 2.1*tolerance:
        return False

    def squared_distance(p):
        foot, radial, _ = _occ_frame(group, p)
        r = np.clip(np.dot(point - foot, radial), group.radius, group.radius + 0.20)
        diff = point - (foot + r*radial)
        return np.dot(diff, diff)

    threshold = tolerance**2
    if any(squared_distance(np.array(on_edge)) <= threshold for on_edge, _ in pairs):
        return True

    def residual(t):
        return squared_distance(np.array(edge.valueAt(t)))

    a, b = edge.ParameterRange
    curve = edge.Curve
    knots = curve.getKnots() if hasattr(curve, "getKnots") else []
    nodes = sorted(set([a, b] + [k for k in knots if a < k < b]))
    for low, high in zip(nodes, nodes[1:]):
        ts = np.linspace(low, high, 13)
        vals = [residual(t) for t in ts]
        if min(vals) <= threshold:
            return True
        for j in range(len(ts)):
            if j and vals[j - 1] < vals[j]:
                continue
            if j + 1 < len(ts) and vals[j + 1] < vals[j]:
                continue
            interval = (ts[max(0, j - 1)], ts[min(len(ts) - 1, j + 1)])
            fit = minimize_scalar(residual, bounds=interval, method="bounded",
                                  options={"xatol": 1e-12})
            if fit.fun <= threshold:
                return True
    return False


def _occ_in_dead_zone(point, stock_boundary, groups, faces, side_edges):
    tol = 0.03
    if stock_boundary.distToShape(Part.Vertex(FreeCAD.Vector(*point)))[0] <= tol:
        return True
    for g, edges in zip(groups, side_edges):
        c, d, r = _occ_frame(g, point)
        if r < g.radius - tol or r > g.radius + 0.20 + tol:
            continue
        if min(abs(r - g.radius), abs(r - g.radius - 0.20)) <= tol:
            if _occ_contains_projection(faces, g, c + g.radius*d):
                return True
        if any(_occ_near_side(point, g, edge, tol) for edge in edges):
            return True
    return False


def verify_occupancy(original, result, groups):
    """Compare majority-ray material membership to the radial projection edit."""
    before, after = _ParitySolid(original), _ParitySolid(result)
    faces = original.Faces
    boundary = Part.makeCompound(original.Shells)
    edge_sets = [_occ_edge_boundary(faces, g) for g in groups]
    tested = 0
    dead_zone_differences = 0
    disagreements = 0
    examples = []
    for point, site in _occ_probes(faces, groups):
        nominal = before.contains(point)
        expected = nominal and not _occ_relief_contains(faces, groups, point)
        observed = after.contains(point)
        tested += 1
        if observed == expected:
            continue
        if _occ_in_dead_zone(point, boundary, groups, faces, edge_sets):
            dead_zone_differences += 1
        else:
            disagreements += 1
            if len(examples) < 8:
                examples.append((site, point.round(7).tolist(), expected, observed))
    print(f"{tested} occupancy probes, {disagreements} disagreements outside "
          f"the 0.03 mm boundary tolerance, {dead_zone_differences} boundary differences",
          flush=True)
    if disagreements:
        raise ValueError(f"Occupancy check failed: {examples}")
    return {"probes": tested, "disagreements": disagreements,
            "boundary_differences": dead_zone_differences}


def mesh_relief_solid(nominal, groups):
    """Regularize intersecting sweeps on a common stock-boundary tessellation."""
    import MeshPart
    import manifold3d as manifold

    # Mesh the entire shell together so adjacent faces share their exact
    # mesh edges. Center large world coordinates before the float mesh API.
    source = nominal.cleaned()
    center = np.asarray(source.BoundBox.Center)
    source.translate(FreeCAD.Vector(*(-center)))
    mesh = MeshPart.meshFromShape(Shape=source, LinearDeflection=0.002,
                                  AngularDeflection=0.15, Relative=False, Segments=True)
    if mesh.countSegments() != len(nominal.Faces):
        raise ValueError("Stock mesh did not retain the original face segments")
    vertices, triangles = mesh.Topology
    points = np.asarray(vertices, dtype=np.float64) + center
    facets = np.asarray(triangles, dtype=np.uint64)
    del vertices, triangles
    def closed_mesh(p, t):
        data = manifold.Mesh64(np.ascontiguousarray(p, dtype=np.float64),
                               np.ascontiguousarray(t, dtype=np.uint64))
        data.merge()
        shape = manifold.Manifold(data)
        if shape.status() != manifold.Error.NoError or shape.is_empty() or shape.volume() <= 0:
            raise ValueError("Could not form a positive closed manifold from the surface triangles")
        return shape

    pieces = [closed_mesh(points, facets)]
    for group in groups:
        for index in group.indices:
            selected = facets[mesh.getSegment(index)]
            local_ids, inverse = np.unique(selected, return_inverse=True)
            local_facets = inverse.reshape(-1, 3)
            near = points[local_ids]
            count = len(near)
            delta = near - group.origin
            radial = delta - np.outer(delta @ group.axis, group.axis)
            far = near + (0.20 / group.radius) * radial

            # Only the boundary half-edges need side walls; an actual periodic
            # seam already has both incident triangles in the source mesh.
            directed = local_facets[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
            _, first, uses = np.unique(np.sort(directed, axis=1), axis=0,
                                       return_index=True, return_counts=True)
            a, b = directed[first[uses == 1]].T
            wall = np.vstack((local_facets, local_facets[:, ::-1] + count,
                              np.column_stack((a, a + count, b + count)),
                              np.column_stack((a, b + count, b))))
            pieces.append(closed_mesh(np.vstack((near, far)), wall))
    result = manifold.Manifold.batch_boolean(pieces, manifold.OpType.Subtract)
    result = result.as_original().simplify(0.002)
    if result.status() != manifold.Error.NoError or sum(c.volume() > 0 for c in result.decompose()) != 1:
        raise ValueError("Regularized relief does not have a single material component")
    triangles = result.to_mesh64()
    point_list = [FreeCAD.Vector(*p) for p in triangles.vert_properties[:, :3]]
    facet_list = [tuple(row) for row in triangles.tri_verts.tolist()]
    shell_shape = Part.Shape()
    shell_shape.makeShapeFromMesh((point_list, facet_list), 1e-7)
    return Part.makeSolid(Part.makeCompound(shell_shape.Shells))


def _one_solid(shape, description):
    solids = shape.Solids
    if len(solids) != 1:
        raise ValueError(f"{description} contains {len(solids)} solids instead of one")
    solid = solids[0]
    if solid.countElement("Face") != shape.countElement("Face"):
        raise ValueError(f"{description} also contains faces outside its solid")
    if not (solid.isClosed() and solid.isValid() and solid.Volume > 0):
        raise ValueError(f"{description} is not a valid, closed, positive-volume solid")
    return solid


def _build(nominal, cutters, sequential):
    # Every tool comes from the original solid. Difference by successive
    # tools and difference by their union describe the same material set.
    result = nominal.copy()
    if sequential:
        for tool in cutters:
            result = result.cut(tool)
    else:
        result = result.cut(tuple(cutters))
    return _one_solid(result, "Cut result")


def issue_fit(part_path, output_dir):
    nominal = _one_solid(Part.read(str(part_path)), "IGES input")
    groups = find_groups(nominal)
    print(f"{part_path.name}: {len(nominal.Faces)} faces, {len(groups)} qualifying groups", flush=True)
    for group in groups:
        numbers = ','.join(str(i + 1) for i in group.indices)
        print(f"  faces {numbers}: r {group.radius:.8f} mm, "
              f"coverage {math.degrees(group.coverage):.8f} deg, span {group.span:.8f} mm", flush=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".relieve-", dir=output_dir) as scratch:
        step = Path(scratch) / "fit.step"

        def publish(candidate):
            candidate.exportStep(str(step))
            reopened = _one_solid(Part.read(str(step)), "Reopened STEP")
            verify_occupancy(nominal, reopened, groups)
            print(f"FreeCAD volume: {nominal.Volume:.8f} -> {reopened.Volume:.8f} mm^3", flush=True)
            os.replace(step, output_dir / "fit.step")

        if not groups:
            publish(nominal)
            return

        try:
            cutters = [tool for group in groups for tool in make_cutters(nominal, group)]
        except (RuntimeError, ValueError, Part.OCCError) as error:
            print(f"Sweep construction failed: {error}", flush=True)
            cutters = []
        if cutters:
            for sequential in ([False, True] if len(cutters) > 1 else [False]):
                try:
                    publish(_build(nominal, cutters, sequential))
                    return
                except (RuntimeError, ValueError, Part.OCCError) as error:
                    print(f"Trimmed-surface Boolean failed validation: {error}", flush=True)

        # Curved cross-drill intersections can defeat the surface-intersection
        # solver. Regularize all the radial sweeps in one common fine mesh,
        # where coincident faces use identical vertices. Export as a B-rep and
        # run the same checks against the original continuous IGES surfaces.
        print("Using a 0.002 mm deflection for the regularized Boolean", flush=True)
        publish(_one_solid(mesh_relief_solid(nominal, groups), "Regularized cut"))


def main():
    parser = argparse.ArgumentParser(description="Increase eligible trimmed IGES bore diameters by 0.40 mm.")
    parser.add_argument("--part", required=True, type=Path, help="Source IGES solid in its original frame")
    parser.add_argument("--out", required=True, type=Path, help="Directory for fit.step")
    args = parser.parse_args()
    issue_fit(args.part, args.out)


if __name__ == "__main__":
    main()
PYTHON_RELIEF_SOURCE
chmod 755 /app/out/relieve.py

/opt/conda/bin/python /app/out/relieve.py --part "/app/inputs/Pump Carriage.iges" --out /app/out

brace_dir=$(mktemp -d)
trap 'rm -rf -- "$brace_dir"' EXIT
/opt/conda/bin/python /app/out/relieve.py --part "/app/inputs/Pump Syringe Brace.iges" --out "$brace_dir"
