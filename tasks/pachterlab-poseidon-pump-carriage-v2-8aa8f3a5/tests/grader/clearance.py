"""Straight insertion against the triangle boundary of the manufactured solid."""

import numpy as np
from numba import njit


@njit
def _hits_disk(polygon, size):
    """Test a convex projected polygon, segment or point against the probe disk."""
    positive = False
    negative = False
    for i in range(size):
        j = (i + 1) % size
        x, y = polygon[i]
        u, v = polygon[j]
        dx, dy = u - x, v - y
        length_squared = dx * dx + dy * dy
        t = 0.0
        if length_squared > 0.0:
            t = min(1.0, max(0.0, -(x * dx + y * dy) / length_squared))
        if (x + t * dx) ** 2 + (y + t * dy) ** 2 <= 1.5 ** 2:
            return True
        side = x * v - y * u
        positive |= side > 0.0
        negative |= side < 0.0
    return size >= 3 and positive != negative


@njit
def _blocked(triangles, low, high, candidates, qx, qy, sign, tip):
    """Clip each potential contact triangle at the tip of the inserted shank."""
    polygon = np.empty((4, 2))
    for index in candidates:
        front = low[index, 2] if sign > 0 else -high[index, 2]
        back = high[index, 2] if sign > 0 else -low[index, 2]
        if back < tip:
            continue
        if front >= tip:
            return True
        tri = triangles[index]
        size = 0
        for i in range(3):
            j = (i + 1) % 3
            a, b = sign * tri[i, 2] - tip, sign * tri[j, 2] - tip
            if a >= 0:
                polygon[size, 0] = tri[i, 0] - qx
                polygon[size, 1] = tri[i, 1] - qy
                size += 1
            if (a > 0 and b < 0) or (a < 0 and b > 0):
                fraction = a / (a - b)
                polygon[size, 0] = tri[i, 0] - qx + fraction * (tri[j, 0] - tri[i, 0])
                polygon[size, 1] = tri[i, 1] - qy + fraction * (tri[j, 1] - tri[i, 1])
                size += 1
        if _hits_disk(polygon, size):
            return True
    return False


@njit
def _approach_limit(triangles, low, high, candidates, qx, qy, sign, first, last):
    """Bisect window stations with the shank extending two millimetres past each."""
    if not _blocked(triangles, low, high, candidates, qx, qy, sign, first - 2.0):
        return first
    if _blocked(triangles, low, high, candidates, qx, qy, sign, last - 2.0):
        return last
    for _ in range(55):
        middle = (first + last) * 0.5
        if middle == first or middle == last:
            break
        if _blocked(triangles, low, high, candidates, qx, qy, sign, middle - 2.0):
            first = middle
        else:
            last = middle
    return (first + last) * 0.5


@njit
def _limits(triangles, low, high, centers, ranges):
    limits = np.empty_like(ranges)
    indices = np.empty(len(triangles), dtype=np.int64)
    polygon = np.empty((3, 2))
    for row, (qx, qy) in enumerate(centers):
        first, last = ranges[row]
        if first == last:
            limits[row, :] = ranges[row, :]
            continue
        count = 0
        for index in range(len(triangles)):
            dx = max(low[index, 0] - qx, qx - high[index, 0], 0.0)
            dy = max(low[index, 1] - qy, qy - high[index, 1], 0.0)
            if dx * dx + dy * dy > 1.5 ** 2:
                continue
            for i in range(3):
                polygon[i, 0] = triangles[index, i, 0] - qx
                polygon[i, 1] = triangles[index, i, 1] - qy
            if _hits_disk(polygon, 3):
                indices[count] = index
                count += 1
        candidates = indices[:count]
        limits[row, 0] = _approach_limit(
            triangles, low, high, candidates, qx, qy, 1, first, last
        )
        limits[row, 1] = -_approach_limit(
            triangles, low, high, candidates, qx, qy, -1, -last, -first
        )
    return limits


def access_rows(triangles, c, d, e, radius, walls):
    """Clip wall intervals to the clear part of each straight probe approach."""
    axes = np.column_stack((e, np.cross(d, e), d))
    local = np.ascontiguousarray((np.asarray(triangles, dtype=float) - c) @ axes)
    angle = np.deg2rad(2.0 * np.arange(180))
    centers = (radius - 1.70) * np.column_stack((np.cos(angle), np.sin(angle)))
    ranges = np.array(
        [(w[0][0], w[-1][1]) if len(w) else (0.0, 0.0) for w in walls], dtype=float
    )
    limits = _limits(local, local.min(axis=1), local.max(axis=1), centers, ranges)
    rows = []
    for wall, (plus_limit, minus_limit) in zip(walls, limits, strict=True):
        plus = [[max(a, plus_limit), b] for a, b in wall if b > plus_limit]
        minus = [[a, min(b, minus_limit)] for a, b in wall if a < minus_limit]
        rows.append([wall, plus, minus])
    return rows
