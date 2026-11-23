"""Tests of the geometric support of indexed triangle boundaries."""

from fractions import Fraction

import numpy as np
import shapely
from trimesh.util import bounds_tree


def _integer_vertices(vertices):
    """Represent the serialized binary coordinates on one exact integer lattice."""
    ratios = [[float(x).as_integer_ratio() for x in row] for row in vertices]
    unit = max(den for row in ratios for _, den in row)
    return [tuple(num * (unit // den) for num, den in row) for row in ratios]


def _normal(a, b, c):
    x, y, z = (j - i for i, j in zip(a, b))
    u, v, w = (j - i for i, j in zip(a, c))
    return y*w - z*v, z*u - x*w, x*v - y*u


def _oriented_area(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _apart(left, right):
    # Both triangles are counterclockwise. Contact on a supporting edge or
    # vertex has empty 2-D interior and is permitted.
    return any(all(_oriented_area(left[i-1], left[i], p) <= 0 for p in right)
               for i in range(3))


def _dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def _plane_interval(triangle, heights, axis):
    """Ends of the triangle/plane intersection along a nonconstant coordinate."""
    points = [p[axis] for p, h in zip(triangle, heights) if h == 0]
    for i in range(3):
        j = i - 1
        if heights[i] * heights[j] < 0:
            a, b = heights[i], heights[j]
            points.append(Fraction(a * triangle[j][axis] - b * triangle[i][axis], a - b))
    return min(points), max(points)


def _face_data(mesh):
    vertices = _integer_vertices(mesh.vertices)
    triangles = [tuple(vertices[i] for i in face) for face in mesh.faces]
    normals = [_normal(*t) for t in triangles]
    heights = [_dot(n, t[0]) for n, t in zip(normals, triangles)]
    boxes = np.column_stack((mesh.triangles.min(axis=1), mesh.triangles.max(axis=1)))
    return triangles, normals, heights, boxes, bounds_tree(boxes)


def coincident_patches(mesh):
    """Find two faces with coplanar area overlap or crossing relative interiors."""
    if not len(mesh.faces):
        return None
    triangles, normals, heights, boxes, tree = _face_data(mesh)
    projections = []
    for normal, triangle in zip(normals, triangles):
        drop = max(range(3), key=lambda i: abs(normal[i]))
        axes = [i for i in range(3) if i != drop]
        tri = tuple(tuple(p[k] for k in axes) for p in triangle)
        projections.append(tri if _oriented_area(*tri) >= 0 else tri[::-1])

    for i, box in enumerate(boxes):
        if not any(normals[i]):
            continue
        for j in tree.intersection(box):
            if j <= i or not any(normals[j]):
                continue
            across_i = [_dot(normals[i], p) - heights[i] for p in triangles[j]]
            if not any(across_i):
                a, b = projections[i], projections[j]
                if not _apart(a, b) and not _apart(b, a):
                    return i, j
                continue
            if not (min(across_i) < 0 < max(across_i)):
                continue
            across_j = [_dot(normals[j], p) - heights[j] for p in triangles[i]]
            if not (min(across_j) < 0 < max(across_j)):
                continue
            # Both triangles straddle the other's plane. They cross strictly
            # inside both faces only if their intervals on the intersection
            # line overlap with positive length.
            direction = _normal((0, 0, 0), normals[i], normals[j])
            axis = max(range(3), key=lambda k: abs(direction[k]))
            a = _plane_interval(triangles[i], across_j, axis)
            b = _plane_interval(triangles[j], across_i, axis)
            if max(a[0], b[0]) < min(a[1], b[1]):
                return i, j
    return None


def buried_face(mesh):
    """Find a triangle with nonzero material winding immediately on its outward side."""
    if not len(mesh.faces):
        return None
    triangles, normals, heights, _, tree = _face_data(mesh)
    # Only the broad phase uses floating centroids. Include summation and
    # division roundoff from the three original vertex coordinates.
    padding = 8 * np.finfo(float).eps * np.abs(mesh.triangles).max(axis=(0, 1))
    for i, (triangle, center, normal) in enumerate(zip(triangles, mesh.triangles_center, normals)):
        k = max(range(3), key=lambda ax: abs(normal[ax]))
        if not normal[k]:
            continue
        direction = 1 if normal[k] > 0 else -1
        u, v = [ax for ax in range(3) if ax != k]
        triple_center = tuple(map(sum, zip(*triangle)))
        box = np.r_[center - padding, center + padding]
        if direction > 0:
            box[3 + k] = mesh.bounds[1, k]
        else:
            box[k] = mesh.bounds[0, k]
        winding = 0
        for j in tree.intersection(box):
            other = normals[j]
            if j == i or not other[k]:
                continue
            sign = 1 if other[k] > 0 else -1
            orientation = -sign if k == 1 else sign
            # Offset the projected point symbolically by (+epsilon,+epsilon²)
            # within the source face. A tangency at a vertex or edge then
            # resolves to an open projected cell without a displacement bar.
            for a, b in zip(triangles[j], triangles[j][1:] + triangles[j][:1]):
                du, dv = b[u] - a[u], b[v] - a[v]
                side = du * (triple_center[v] - 3*a[v]) - dv * (triple_center[u] - 3*a[u])
                if (side or -dv or du) * orientation < 0:
                    break
            else:
                # Compare exact depths along the shifted ray. Its origin moves
                # on the source plane and its positive direction is outward.
                depth = (3 * heights[j] - _dot(other, triple_center)) * normal[k]
                slope_u = normal[u] * other[k] - other[u] * normal[k]
                slope_v = normal[v] * other[k] - other[v] * normal[k]
                if (depth or slope_u or slope_v) * sign > 0:
                    winding += direction * sign
        if winding:
            return f'Face {i} has oriented winding {winding} on its outward side'
    return None


def _quantization(mesh):
    """Half a coordinate ULP at the largest magnitude encoded on each axis."""
    scale = np.abs(mesh.triangles).max(axis=(0, 1))
    types = mesh.metadata.get('coordinate_dtypes', (np.float64,) * 3)
    return np.array([abs(float(np.spacing(np.array(x, dtype=t)))) / 2
                     if np.issubdtype(np.dtype(t), np.floating) else 0.
                     for x, t in zip(scale, types)])


def _decimal_steps(mesh):
    """Local ASCII rounding uncertainty constrained by the coordinate column's precision."""
    steps = mesh.metadata.get('ascii_coordinate_half_steps')
    if steps is None:
        return np.zeros_like(mesh.vertices)
    radii = np.asarray(steps).copy()
    coordinates = np.abs(mesh.vertices)
    # A short token such as "33.1" can be a trimmed value in a column that
    # otherwise carries six significant digits or six fixed decimal places.
    # It must not give the whole facet a 0.05 mm uncertainty.
    for axis in range(3):
        radius = radii[:, axis]
        size = coordinates[:, axis]
        nonzero = (size > 0) & (radius > 0)
        if not nonzero.any():
            continue
        decimal_place = np.rint(np.log10(2 * radius[nonzero]))
        magnitude = np.floor(np.log10(size[nonzero]))
        digits = int(np.max(magnitude - decimal_place + 1))
        fixed = radius[radius > 0].min()
        significant = np.zeros_like(radius)
        positive = size > 0
        significant[positive] = 0.5 * 10.**(np.floor(np.log10(size[positive])) + 1 - digits)
        radii[:, axis] = np.minimum(radius, np.maximum(fixed, significant))
    return radii


def _matching_triangles(left, right):
    # Reindexing, cyclic order and element order change no point of a triangle.
    def canonical(triangles):
        return sorted(tuple(sorted(map(tuple, triangle))) for triangle in triangles)
    return len(left.faces) == len(right.faces) and canonical(left.triangles) == canonical(right.triangles)


def _clip(points, heights, top):
    """Clip a convex polygon to an affine half-space with height <= top."""
    within = heights <= top
    if within.all():
        return points, heights
    clipped, levels = [], []
    for i in range(len(heights)):
        j = i - 1
        if within[i] != within[j]:
            t = (top - heights[j]) / (heights[i] - heights[j])
            clipped.append(points[j] + t * (points[i] - points[j]))
            levels.append(top)
        if within[i]:
            clipped.append(points[i])
            levels.append(heights[i])
    return np.asarray(clipped).reshape(-1, 2), np.asarray(levels)


def _uncovered(source, target, precision, source_steps, target_steps):
    """Find a source triangle with points outside the target's serialization envelope."""
    targets = target.triangles
    target_radii = target_steps[target.faces].max(axis=1)
    tree = bounds_tree(np.stack((targets.min(axis=1) - target_radii,
                                 targets.max(axis=1) + target_radii), axis=1))
    exact = _integer_vertices(source.vertices)
    for face_id, (indices, triangle) in enumerate(zip(source.faces, source.triangles)):
        # The dominant-axis projection is nonsingular for any nonzero-area
        # triangle. The exact normal also works for extremely slender facets.
        normal = _normal(*(exact[i] for i in indices))
        drop = max(range(3), key=lambda k: abs(normal[k]))
        if normal[drop] == 0:
            continue
        axes = [i for i in range(3) if i != drop]
        ratios = np.array([normal[k] / normal[drop] for k in axes])
        radius = precision + source_steps[indices].max(axis=0)
        box = np.r_[triangle.min(axis=0) - radius, triangle.max(axis=0) + radius]
        ids = list(tree.intersection(box))
        if not ids:
            return face_id
        local = targets[ids] - triangle[0]
        xy = local[:, :, axes]
        heights = local[:, :, drop] + xy @ ratios
        allowance = radius + target_radii[ids]
        band = allowance[:, drop] + allowance[:, axes] @ abs(ratios)
        candidates = (heights.min(axis=1) <= band) & (heights.max(axis=1) >= -band)
        corners = np.array([[-1., -1.], [-1., 1.], [1., 1.], [1., -1.]])
        patches = []
        for j in np.flatnonzero(candidates):
            points, h = _clip(xy[j], heights[j], band[j])
            points, _ = _clip(points, -h, band[j])
            shifts = corners * allowance[j, axes]
            if len(points):
                expanded = (points[:, None, :] + shifts).reshape(-1, 2)
                patches.append(shapely.convex_hull(shapely.multipoints(expanded)))
        footprint = shapely.Polygon((triangle - triangle[0])[:, axes])
        remainder = footprint
        for patch in sorted(patches, key=lambda p: -p.area):
            remainder = shapely.difference(remainder, patch)
            if remainder.is_empty:
                break
        if not remainder.is_empty:
            return face_id
    return None


def surface_difference(left, right):
    """Return a witness if either complete triangle surface lacks coverage by the other."""
    if not len(left.faces) or not len(right.faces):
        return 'one study has an empty surface' if len(left.faces) != len(right.faces) else None
    if _matching_triangles(left, right):
        return None
    # Independent coordinate rounding in the two encodings adds. The projected
    # slab and footprint calculations use local differences, two products, a
    # dot sum, and interpolation. gamma(16) bounds their double roundoff at
    # the coordinate scale, independently of the production dimension budget.
    scale = np.maximum(np.abs(left.triangles).max(axis=(0, 1)), np.abs(right.triangles).max(axis=(0, 1)))
    eps = np.finfo(np.float64).eps / 2
    precision = _quantization(left) + _quantization(right) + scale * (16 * eps / (1 - 16 * eps))
    left_steps, right_steps = _decimal_steps(left), _decimal_steps(right)
    for name, source, target, sr, tr in [('first', left, right, left_steps, right_steps),
                                         ('second', right, left, right_steps, left_steps)]:
        face = _uncovered(source, target, precision, sr, tr)
        if face is not None:
            return f'{name} study face {face} has surface outside serialization precision of its counterpart'
    return None
