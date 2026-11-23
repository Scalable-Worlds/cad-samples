"""Compare material at signed surface, radial and trim stations of the original stock."""
import io
import math
import os
from collections import defaultdict

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import FreeCAD
import Part
import numpy as np
from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Builder
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepTools import BRepTools
from OCP.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCP.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCP.IntCurveSurface import IntCurveSurface_Tangent
from OCP.Precision import Precision
from OCP.TopAbs import TopAbs_FACE, TopAbs_IN
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Compound, TopoDS_Shape
from OCP.gp import gp_Dir, gp_Lin, gp_Pnt, gp_Pnt2d
from scipy.spatial import cKDTree
from scipy.stats import qmc

from supports import LIMIT_EPS

T = 0.03
ALLOWANCE = 0.20
RAYS = [(0.33585, 0.7985, 0.50049), (0.59211, -0.18258, 0.78342),
        (-0.56118, 0.61272, 0.56143)]
CATEGORIES = ("normal_stock", "normal_void", "bore_void", "allowance", "backing",
              "axial_trim", "hoop_trim", "partial_support", "radius_excluded",
              "export_boundary", "nominal_material", "span_excluded")


def xyz(vector):
    return np.asarray(tuple(vector), dtype=float)


def radial_system(point, origin, axis):
    delta = point - origin
    center = origin + np.dot(delta, axis) * axis
    radial = point - center
    return center, radial / np.linalg.norm(radial)


def stations(faces, groups):
    # Both requested sample sizes use the same prefix of the unscrambled sequence.
    unit_uv = qmc.Halton(2, scramble=False).random(32 * 16 + 2)[1:]
    by_face = {i: group for group in groups if group.eligible for i in group.faces}
    sites, labels, face_numbers = [], [], []

    def add(point, label, index):
        sites.append(point)
        labels.append(label)
        face_numbers.append(index)

    for index, face in enumerate(faces):
        u0, u1, v0, v1 = face.ParameterRange
        trim = BRepTopAdaptor_FClass2d(TopoDS.Face(occ_shape(face)), Precision.PConfusion_s())
        g = by_face.get(index)
        if g is not None:
            origin, axis = xyz(g.origin), xyz(g.axis)
        n = 16 if g is not None else 8
        count = 0
        for seq_index, (a, b) in enumerate(unit_uv[:32 * n + 1]):
            u, v = u0 + a * (u1 - u0), v0 + b * (v1 - v0)
            if trim.Perform(gp_Pnt2d(u, v)) != TopAbs_IN:
                continue
            point = xyz(face.valueAt(u, v))
            if count < 8 and seq_index < 32 * 8 + 1:
                normal = xyz(face.normalAt(u, v))
                for offset in (-0.07, 0.07):
                    add(point + offset * normal, "normal", index)
            if g is not None:
                center, radial = radial_system(point, origin, axis)
                for offset in (-0.07, 0.06, 0.12, 0.16, ALLOWANCE - 1.08 * T,
                               ALLOWANCE + 1.08 * T, 0.24, 0.28):
                    kind = "bore_void" if offset < 0 else "backing" if offset > ALLOWANCE else "allowance"
                    add(center + (g.radius + offset) * radial, kind, index)
            count += 1
            if count == n:
                break
        if g is None:
            continue
        for edge in face.Edges:
            if edge.isSeam(face):
                continue
            lo, hi = edge.ParameterRange
            for fraction in (0.15, 0.50, 0.85):
                point = xyz(edge.valueAt(lo + fraction * (hi - lo)))
                center, radial = radial_system(point, origin, axis)
                hoop = np.cross(axis, radial)
                for distance in (-0.12, -1.5 * T, 1.5 * T, 0.12):
                    angle = distance / g.radius
                    turned = math.cos(angle) * radial + math.sin(angle) * hoop
                    for offset in ((0.10, 0.25) if abs(distance) == 0.12 else (0.10,)):
                        radius = g.radius + offset
                        add(center + distance * axis + radius * radial, "axial_trim", index)
                        add(center + radius * turned, "hoop_trim", index)
    return np.asarray(sites), np.array(labels), np.array(face_numbers)


def occ_shape(shape):
    result = TopoDS_Shape()
    BRepTools.Read_s(result, io.BytesIO(shape.exportBrepToString().encode()), BRep_Builder())
    return result


def face_batches(shape):
    # Partition faces spatially. Intersections are still OCCT shape-intersector
    # crossings, with coincident hits merged across all batches on the same ray.
    builder = BRep_Builder()
    faces, boxes = [], []
    walk = TopExp_Explorer(occ_shape(shape), TopAbs_FACE)
    while walk.More():
        face = walk.Current()
        box = Bnd_Box()
        BRepBndLib.Add_s(face, box, False)
        faces.append(face)
        boxes.append(box.Get())
        walk.Next()
    boxes = np.asarray(boxes)
    centers = (boxes[:, :3] + boxes[:, 3:]) * 0.5
    pending = [np.arange(len(faces))]
    while pending:
        indices = pending.pop()
        if len(indices) > 8:
            axis = np.argmax(np.ptp(centers[indices], axis=0))
            ordered = indices[np.argsort(centers[indices, axis])]
            middle = len(ordered) // 2
            pending.extend([ordered[:middle], ordered[middle:]])
            continue
        compound = TopoDS_Compound()
        builder.MakeCompound(compound)
        for i in indices:
            builder.Add(compound, faces[i])
        local = boxes[indices]
        yield compound, local[:, :3].min(axis=0)-1e-7, local[:, 3:].max(axis=0)+1e-7


def inside(shape, sites):
    # Two matching parities already decide the three-ray majority.
    batches = list(face_batches(shape))
    votes = np.zeros((len(RAYS), len(sites)), dtype=bool)
    active = np.arange(len(sites))
    casters = []
    for compound, lo, hi in batches:
        it = IntCurvesFace_ShapeIntersector()
        it.Load(compound, 1e-7)
        casters.append((it, lo, hi))
    for ray_index, dr in enumerate(RAYS):
        if ray_index == 2:
            active = np.where(votes[0] != votes[1])[0]
            if not len(active):
                break
        pts = sites[active]
        direction = np.array(dr) / np.linalg.norm(dr)
        gp_direction = gp_Dir(*dr)
        hits = [[] for _ in pts]
        lines = [gp_Lin(gp_Pnt(*p), gp_direction) for p in pts]
        for it, lo, hi in casters:
            a, b = (lo - pts) / direction, (hi - pts) / direction
            entry = np.maximum(0.0, np.minimum(a, b).max(axis=1))
            leave = np.maximum(a, b).min(axis=1)
            for j in np.flatnonzero(leave >= entry):
                it.Perform(lines[j], 0.0, 1.e100)
                assert it.IsDone(), "OCCT could not compute a material crossing"
                for k in range(1, it.NbPnt() + 1):
                    distance = it.WParameter(k)
                    if distance > 0 and it.Transition(k) != IntCurveSurface_Tangent:
                        hits[j].append(distance)
        for j, crossings in zip(active, hits):
            sorted_hits = sorted(crossings)
            distinct = []
            for w in sorted_hits:
                if not distinct or w - distinct[-1] > 1e-6:
                    distinct.append(w)
            votes[ray_index, j] = len(distinct) % 2 == 1
    return (votes[0] & votes[1]) | ((votes[0] | votes[1]) & votes[2])


def aabb(shape):
    # A triangulation box can miss a cylindrical extremum between vertices.
    box = Bnd_Box()
    BRepBndLib.Add_s(occ_shape(shape), box, False)
    bounds = np.asarray(box.Get())
    return bounds[:3], bounds[3:]


def nearby(point, bounds, margin):
    lo, hi = bounds
    return np.all(point >= lo-margin) and np.all(point <= hi+margin)


class Band:
    def __init__(self, group, faces):
        self.origin = xyz(group.origin)
        self.axis = xyz(group.axis)
        self.radius = group.radius
        self.faces = [faces[i] for i in group.faces]
        self.bounds = [aabb(f) for f in self.faces]
        self.outer_and_sides = None

    def projection(self, point):
        delta = point - self.origin
        center = self.origin + np.dot(delta, self.axis) * self.axis
        ray = point - center
        rho = np.linalg.norm(ray)
        if rho == 0:
            return rho, None
        return rho, center + (self.radius/rho) * ray

    def contains(self, point):
        rho, q = self.projection(point)
        if not self.radius < rho < self.radius + ALLOWANCE:
            return False
        # Measure to the trimmed face itself. This handles holes, periodic
        # carriers and UV projections that would otherwise clamp to a surface end.
        vertex = None
        for face, bounds in zip(self.faces, self.bounds):
            if nearby(q, bounds, 1e-5):
                if vertex is None:
                    vertex = Part.Vertex(*q)
                if face.distToShape(vertex)[0] <= 1e-5:
                    return True
        return False

    def make_outer_and_sides(self):
        # Radial scaling about a line maps an entire cylindrical patch exactly
        # to its outer cap. The fibers over its true boundary edges are ruled.
        affine = np.eye(4)
        dilation = (ALLOWANCE / self.radius) * (np.eye(3) - np.outer(self.axis, self.axis))
        affine[:3,:3] += dilation
        affine[:3,3] = -dilation @ self.origin
        matrix = FreeCAD.Matrix(*affine.ravel())
        boundaries = []
        edge_buckets = defaultdict(list)
        for face in self.faces:
            boundaries.append(face.transformGeometry(matrix))
            for edge in face.Edges:
                if edge.isSeam(face):
                    continue
                edge_buckets[edge.hashCode()].append(edge)
        for bucket in edge_buckets.values():
            for edge in bucket:
                if sum(edge.isSame(other) for other in bucket) == 1:
                    outer = edge.transformGeometry(matrix)
                    boundaries.append(Part.makeRuledSurface(edge, outer))
        self.outer_and_sides = [(item, aabb(item)) for item in boundaries]

    def boundary_distance(self, point, vertex, search):
        rho, _ = self.projection(point)
        if not self.radius - search <= rho <= self.radius + ALLOWANCE + search:
            return math.inf
        if self.outer_and_sides is None:
            self.make_outer_and_sides()
        # The inner cap is already among the original stock faces.
        return min((face.distToShape(vertex)[0] for face, bounds in self.outer_and_sides
                    if nearby(point, bounds, search)), default=math.inf)


class BoundaryDistances:
    def __init__(self, faces, bands, sample_points):
        self.faces = faces
        self.bands = bands
        self.bounds = np.array([aabb(face) for face in faces])
        self.samples = cKDTree(sample_points)

    def local_distance(self, point, search):
        # A distance above the search radius can be represented by infinity.
        # Inside the tolerance zone the exact distance is immaterial.
        sample_distance = self.samples.query(point, distance_upper_bound=T)[0]
        if sample_distance <= T:
            return sample_distance
        vertex = Part.Vertex(*point)
        intersects = ((point >= self.bounds[:, 0] - search)
                      & (point <= self.bounds[:, 1] + search))
        best = math.inf
        for j in np.flatnonzero(np.all(intersects, axis=1)):
            best = min(best, self.faces[j].distToShape(vertex)[0])
            if best <= T:
                return best
        for band in self.bands:
            best = min(best, band.boundary_distance(point, vertex, search))
            if best <= T:
                return best
        return best

    def protected_step(self, point):
        clearance = self.local_distance(point, 1.5 * T)
        return min(T / 4, max(0.0, (clearance - T) / 2))


def face_samples(face, count):
    u0, u1, v0, v1 = face.ParameterRange
    trim = BRepTopAdaptor_FClass2d(TopoDS.Face(occ_shape(face)), Precision.PConfusion_s())
    sequence = qmc.Halton(2, scramble=False).random(32 * count + 2)[1:]
    seen = 0
    for a, b in sequence:
        u, v = u0 + a * (u1 - u0), v0 + b * (v1 - v0)
        if trim.Perform(gp_Pnt2d(u, v)) == TopAbs_IN:
            yield xyz(face.valueAt(u, v)), xyz(face.normalAt(u, v))
            seen += 1
            if seen == count:
                break


def export_samples(shape, original_face_count, boundary):
    # Sample every face of ordinary exchanges. For a densely faceted exchange,
    # stratify by surface area to bound ray work independently of tessellation.
    faces = shape.Faces
    budget = 4 * original_face_count
    if len(faces) <= budget:
        indices = range(len(faces))
    else:
        areas = np.cumsum([f.Area for f in faces])
        targets = (np.arange(budget) + 0.5) * (areas[-1] / budget)
        indices = np.unique(np.searchsorted(areas, targets))
    for i in indices:
        for point, normal in face_samples(faces[i], 2):
            yield point - 0.07 * normal, i
            yield point + 0.07 * normal, i
            # Where the output boundary penetrates protected material or void,
            # sample *within* that penetration on either side of the output
            # face. Fixed offsets alone can jump across a narrow discrepancy.
            distance = boundary.protected_step(point)
            if distance:
                yield point - distance * normal, i
                yield point + distance * normal, i


def material_difference_samples(original, fitted, boundary):
    # Compare the complete B-reps in both directions. The Boolean solids locate
    # changed regions, even if buried. A change fails only on a parity witness
    # outside the original surface's tolerance zone. An approximated skin, for
    # instance, can have many nonempty Boolean fragments and still be compliant.
    for left, right in ((original, fitted), (fitted, original)):
        try:
            difference = left.cut(right)
        except Part.OCCError:
            # Coincident and faceted surfaces can defeat the Boolean. Probe the
            # left operand directly on both sides of its boundary in this direction.
            for point, _ in export_samples(left, len(original.Faces), boundary):
                yield point
            continue
        components = difference.Solids
        if not components:
            # Identical regions give an empty difference. Confirm parity inside
            # the operand's bounds without requiring any nonzero missing volume.
            yield xyz(left.CenterOfMass)
        for component in components:
            yield xyz(component.CenterOfMass)
            for face in component.Faces:
                for point, normal in face_samples(face, 2):
                    yield point - 0.07 * normal
                    distance = boundary.protected_step(point)
                    if distance:
                        yield point - distance * normal


def audit(shape, groups, fitted):
    """Compute the independent material verdict at every published station."""
    faces = shape.Faces
    sites, tags, source_indices = stations(faces, groups)
    bands = [Band(g, faces) for g in groups if g.eligible]
    # Midpoints of opposite normal probes lie exactly on the nominal faces.
    # Index them and their corresponding outer-cap samples before requesting
    # an exact point/face distance for an exchanged surface.
    normal = tags == "normal"
    marks = sites[normal].reshape(-1, 2, 3).mean(axis=1)
    owners = source_indices[normal][::2]
    samples = [marks]
    for group in groups:
        if group.eligible:
            p = marks[np.isin(owners, group.faces)]
            delta = p - group.origin
            radial = delta - np.outer(delta @ group.axis, group.axis)
            samples.append(p + (ALLOWANCE / group.radius) * radial)
    boundary = BoundaryDistances(faces, bands, np.vstack(samples))
    extra = [(p, "export_boundary", i) for p, i in export_samples(fitted, len(faces), boundary)]
    if not bands:
        extra.extend((p, "nominal_material", -1) for p in material_difference_samples(shape, fitted, boundary))
    if extra:
        positions, labels, indices = zip(*extra)
        sites = np.vstack((sites, positions))
        tags = np.concatenate((tags, labels))
        source_indices = np.concatenate((source_indices, indices))
    points, inverse = np.unique(sites, axis=0, return_inverse=True)
    original = inside(shape, points)
    expected = original.copy()
    for band in bands:
        radii = np.linalg.norm(np.cross(points - band.origin, band.axis), axis=1)
        near_ring = ((radii >= band.radius - LIMIT_EPS)
                     & (radii <= band.radius + ALLOWANCE + LIMIT_EPS))
        for i in np.flatnonzero(expected & near_ring):
            if band.contains(points[i]):
                expected[i] = False
    original = original[inverse]
    expected = expected[inverse]
    actual = inside(fitted, points)[inverse]
    preserved = {}
    for group in groups:
        radius_in_range = 1.25 - LIMIT_EPS <= group.radius <= 10 + LIMIT_EPS
        enough_azimuth = group.coverage >= math.pi - LIMIT_EPS
        if group.span < 1 - LIMIT_EPS:
            if radius_in_range and enough_azimuth:
                preserved.update(dict.fromkeys(group.faces, "span_excluded"))
            continue
        if radius_in_range and not enough_azimuth:
            preserved.update(dict.fromkeys(group.faces, "partial_support"))
        if enough_azimuth and not radius_in_range:
            preserved.update(dict.fromkeys(group.faces, "radius_excluded"))

    tags = list(tags)
    for i, tag in enumerate(tags):
        if tag == "normal":
            if expected[i] and source_indices[i] in preserved:
                tags[i] = preserved[source_indices[i]]
            else:
                tags[i] = "normal_stock" if original[i] else "normal_void"
    report = {key: {"tested": 0, "errors": 0, "first": None} for key in CATEGORIES}
    for tag in tags:
        report[tag]["tested"] += 1
    for i in np.flatnonzero(actual != expected):
        point = sites[i]
        if boundary.local_distance(point, T) <= T:
            continue
        stats = report[tags[i]]
        stats["errors"] += 1
        if stats["first"] is None:
            desc = "material" if expected[i] else "void"
            source = f"face {source_indices[i]+1}" if tags[i] != "nominal_material" else "whole solid"
            stats["first"] = f"{source}, point {point.round(5).tolist()}, expected {desc}"
    return report
