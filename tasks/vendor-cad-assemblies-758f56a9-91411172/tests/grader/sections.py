"""Line sections through triangle boundaries in orthonormal pull coordinates."""

import numpy as np
import shapely


CROSSING_MERGE = 1e-6
EMPTY = np.empty((0, 2), dtype=float)


def _without_segment(intervals, low, high):
    """Remove a run along the boundary from a solid-interior section."""
    pieces = []
    for a, b in intervals:
        if a < min(b, low):
            pieces.append((a, min(b, low)))
        if max(a, high) < b:
            pieces.append((max(a, high), b))
    return np.asarray(pieces).reshape(-1, 2)


def _parallel_span(triangle, point, precision):
    """A pull-parallel boundary triangle intersects the ray in a segment."""
    uv = triangle[:, :2]
    axis = int(np.argmax(np.ptp(uv, axis=0)))
    s = uv[:, axis]
    a, b = np.argmin(s), np.argmax(s)
    if s[b] == s[a]:
        return None
    along = uv[b] - uv[a]
    delta = point - uv[a]
    if abs(along[0] * delta[1] - along[1] * delta[0]) > precision * np.abs(along).sum():
        return None
    target = point[axis]
    heights = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        if target < min(s[i], s[j]) - precision or target > max(s[i], s[j]) + precision:
            continue
        if abs(s[i] - s[j]) <= precision:
            heights.extend(triangle[[i, j], 2])
        else:
            t = np.clip((target - s[i]) / (s[j] - s[i]), 0.0, 1.0)
            heights.append(triangle[i, 2] + t * (triangle[j, 2] - triangle[i, 2]))
    return (min(heights), max(heights)) if heights else None


def _winding_sections(triangles, point, heights):
    """Resolve an ambiguous edge hit by testing the open runs between events."""
    levels = np.sort(heights)
    if len(levels) < 2:
        return EMPTY
    groups = np.r_[0, np.flatnonzero(np.diff(levels) > CROSSING_MERGE) + 1]
    levels = np.add.reduceat(levels, groups) / np.diff(np.r_[groups, len(levels)])
    regions = []
    for low, high in zip(levels[:-1], levels[1:]):
        origin = (*point, 0.5 * (low + high))
        r = triangles - origin
        # Signed solid angle of each triangle seen from the sample.
        r /= np.maximum(np.linalg.norm(r, axis=2, keepdims=True), np.finfo(float).tiny)
        a, b, c = r[:, 0], r[:, 1], r[:, 2]
        numerator = np.einsum("ij,ij->i", a, np.cross(b, c))
        denominator = 1 + np.einsum("ij,ij->i", a, b) + np.einsum("ij,ij->i", b, c) + np.einsum("ij,ij->i", c, a)
        if abs(np.arctan2(numerator, denominator).sum()) > np.pi:
            if regions and regions[-1][1] == low:
                regions[-1][1] = high
            else:
                regions.append([low, high])
    return np.asarray(regions).reshape(-1, 2)


def solid_sections(projected_triangles, points, *, bounded_air=False):
    """Return open sections of the solid, or of the air trapped in its pull envelope."""
    triangles = np.asarray(projected_triangles, dtype=float).reshape(-1, 3, 3)
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    result = [EMPTY] * len(points)
    if not len(triangles) or not len(points):
        return result

    # Use original orientations only for the fallback winding-number test.
    triangles_ccw = triangles.copy()
    uv = triangles_ccw[:, :, :2]
    e, f = uv[:, 1] - uv[:, 0], uv[:, 2] - uv[:, 0]
    area2 = e[:, 0] * f[:, 1] - e[:, 1] * f[:, 0]
    reverse = area2 < 0
    triangles_ccw[reverse] = triangles_ccw[reverse][:, [0, 2, 1]]
    precision = 16 * np.finfo(float).eps * max(1., np.max(np.abs(uv)))
    area2 = np.abs(area2)
    flat = area2 <= precision * np.ptp(uv, axis=1).sum(axis=1)
    uv1 = np.roll(uv, -1, axis=1)
    edge = uv1 - uv
    # Canonical arithmetic along a shared edge gives adjacent triangles the
    # same zero predicate. Its half-open rule uses a symbolic shift (+eps,+eps²).
    forward = (edge[:, :, 0] > 0) | ((edge[:, :, 0] == 0) & (edge[:, :, 1] >= 0))
    start = np.where(forward[:, :, None], uv, uv1)
    sign = np.where(forward, 1., -1.)
    direction = edge * sign[:, :, None]
    owns = (edge[:, :, 1] < 0) | ((edge[:, :, 1] == 0) & (edge[:, :, 0] > 0))
    edge_error = precision * np.abs(edge).sum(axis=2)
    low, high = uv.min(axis=1) - precision, uv.max(axis=1) + precision
    tree = shapely.STRtree(shapely.box(low[:, 0], low[:, 1], high[:, 0], high[:, 1]))
    block = max(1, min(2048, 2000000 // len(triangles)))
    for first in range(0, len(points), block):
        q = points[first:first + block]
        i, j = tree.query(shapely.points(q))
        horizontal = flat[j]
        spans = {}
        for p, face in zip(i[horizontal], j[horizontal]):
            span = _parallel_span(triangles[face], q[p], precision)
            if span is not None and span[1] > span[0]:
                spans.setdefault(p, []).append(span)
        p, face = i[~horizontal], j[~horizontal]
        delta = q[p, None, :] - start[face]
        dxy = direction[face]
        side = sign[face] * (dxy[:, :, 0] * delta[:, :, 1] - dxy[:, :, 1] * delta[:, :, 0])
        near = np.abs(side) <= edge_error[face]
        inclusive = np.all((side > 0) | near, axis=1)
        p, face, side, near = p[inclusive], face[inclusive], side[inclusive], near[inclusive]
        # Weights opposite the three directed edges correspond to corners 2,0,1.
        side[near] = 0.
        normalizer = side.sum(axis=1)
        regular = normalizer > 0
        weights = side[regular] / normalizer[regular, None]
        z = (weights * triangles_ccw[face[regular]][:, [2, 0, 1], 2]).sum(axis=1)
        p, owned = p[regular], np.all(~near[regular] | owns[face[regular]], axis=1)
        ray, height = p[owned], z[owned]
        odd_rows = []
        if len(ray):
            order = np.lexsort((height, ray))
            ray, height = ray[order], height[order]
            groups = np.r_[0, np.flatnonzero((ray[1:] != ray[:-1]) | (np.diff(height) > CROSSING_MERGE)) + 1]
            count = np.diff(np.r_[groups, len(height)])
            odd = (count % 2) == 1
            ray = ray[groups][odd]
            height = (np.add.reduceat(height, groups) / count)[odd]
            # Each crossing flips parity. An even multiplicity is a tangency.
            if len(ray):
                rows, begin = np.unique(ray, return_index=True)
                end = np.r_[begin[1:], len(ray)]
                for row, a, b in zip(rows, begin, end):
                    if (b - a) % 2:
                        odd_rows.append(row)
                    else:
                        result[first + row] = height[a:b].reshape(-1, 2)
        for row in odd_rows:
            levels = list(z[p == row])
            for span in spans.get(row, []):
                levels.extend(span)
            result[first + row] = _winding_sections(triangles, q[row], levels)
        for row, cuts in spans.items():
            for a, b in cuts:
                result[first + row] = _without_segment(result[first + row], a, b)
        if bounded_air:
            last = first + len(q)
            result[first:last] = trapped_sections(result[first:last])
            # A pull-parallel face can lie between two interior stock runs.
            # That face belongs to the stock boundary, outside the open steel.
            for row, cuts in spans.items():
                for a, b in cuts:
                    result[first + row] = _without_segment(result[first + row], a, b)
    return result


def trapped_sections(stock):
    """Return exactly the bounded air between consecutive stock intervals."""
    return [
        np.column_stack((segments[:-1, 1], segments[1:, 0])) if len(segments) > 1 else EMPTY
        for segments in stock
    ]


def difference_lengths(left, right):
    """Lebesgue lengths of the symmetric differences of two interval lists."""
    values = np.zeros(len(left))
    for row, (a, b) in enumerate(zip(left, right, strict=True)):
        if not len(a) and not len(b):
            continue
        total = (a[:, 1] - a[:, 0]).sum() + (b[:, 1] - b[:, 0]).sum()
        i = j = 0
        while i < len(a) and j < len(b):
            a0, a1 = a[i]
            b0, b1 = b[j]
            total -= 2 * max(0., min(a1, b1) - max(a0, b0))
            if a1 <= b1:
                i += 1
            else:
                j += 1
        values[row] = max(0., total)
    return values
