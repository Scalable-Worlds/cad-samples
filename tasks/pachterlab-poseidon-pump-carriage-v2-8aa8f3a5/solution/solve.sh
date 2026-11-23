#!/bin/bash
set -euo pipefail

out=/app/out
mkdir -p "$out"
rm -f -- "$out/bore_access.py" "$out/base/survey.json" "$out/screen-leg/survey.json" "$out/top-stage/survey.json"
cat > "$out/bore_access.py" <<'BORE_ACCESS_PY'
#!/opt/conda/bin/python
"""Compute trimmed cylindrical wall and straight-insertion probe surveys from a binary STL."""

import argparse
import json
import math
import os
import struct
import tempfile
from pathlib import Path

import numpy as np
from numba import njit
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def read_triangles(path):
    """Read the triangles in their stored winding (binary STL, millimetres)."""
    record = np.dtype([('normal','<f4',3), ('vertices','<f4',(3,3)), ('attribute','<u2')])
    with open(path, 'rb') as stream:
        header = stream.read(84)
        if len(header) != 84:
            raise ValueError('incomplete STL header')
        count, = struct.unpack_from('<I', header, 80)
        if count == 0 or os.fstat(stream.fileno()).st_size != 84 + record.itemsize * count:
            raise ValueError('invalid binary STL triangle count or file size')
        facets = np.fromfile(stream, dtype=record, count=count)
    if len(facets) != count or not np.isfinite(facets['vertices']).all():
        raise ValueError('invalid STL vertex data')
    return np.ascontiguousarray(facets['vertices'], dtype=np.float64)


def _cyl_dot(x, y):
    return np.einsum('...i,...i->...', x, y)


def _cyl_frame(d):
    u = np.cross(d, np.eye(3)[np.argmin(np.abs(d))])
    u /= np.linalg.norm(u)
    return np.array([u, np.cross(d, u)])


def _land_adjacency(faces, vertex_count):
    edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    keys = edges[:, 0] * vertex_count + edges[:, 1]
    order = np.argsort(keys)
    pairs = keys[order[:-1]] == keys[order[1:]]
    return edges, order[:-1][pairs], order[1:][pairs]


def _crease_circles(vertices, faces, normals, adjacency):
    """A concave change of chord plane gives a generator and three circle points."""
    edges, first, second = adjacency
    i, j = first // 3, second // 3
    opposite1 = faces[i, (first + 2) % 3]
    opposite2 = faces[j, (second + 2) % 3]
    a, b = vertices[edges[first]].transpose(1, 0, 2)
    x, y = vertices[opposite1] - a, vertices[opposite2] - a
    cosine = _cyl_dot(normals[i], normals[j])
    h = np.linalg.norm(b - a, axis=1)
    use = ((cosine < 1 - 1e-10) & (cosine > np.cos(np.deg2rad(18.))) &
           (_cyl_dot(y, normals[i]) > 1e-8) & (h > 1e-6))
    a, b, x, y, h = (z[use] for z in (a, b, x, y, h))
    d = (b - a) / h[:, None]
    x -= _cyl_dot(x, d)[:, None] * d
    y -= _cyl_dot(y, d)[:, None] * d
    width = np.linalg.norm(x, axis=1)
    u = x / np.maximum(width, 1e-30)[:, None]
    w = np.cross(d, u)
    denom = 2 * _cyl_dot(w, y)
    good = np.abs(denom) > 1e-10
    offset = np.zeros_like(denom)
    offset[good] = (_cyl_dot(y, y)[good] - width[good] * _cyl_dot(y, u)[good]) / denom[good]
    center = a + 0.5 * x + offset[:, None] * w
    radius = np.hypot(0.5 * width, offset)
    keep = good & (radius >= 1.65) & (radius <= 20.1)
    # Center on the axis near the (already centered) mesh.
    center -= _cyl_dot(center, d)[:, None] * d
    return d[keep], center[keep], radius[keep], h[keep]


def _crease_seeds(d, center, radius, height):
    if not len(radius):
        return
    # dd^T is identical for the two signs of the same axis. This avoids
    # dividing a cluster at an arbitrary hemisphere boundary after rotation.
    i, j = np.triu_indices(3)
    keys = np.column_stack((40 * d[:, i] * d[:, j], center, radius))
    tree = cKDTree(keys)
    # Start fine enough to preserve nearby nominal radii. Also gather wider
    # clusters for the noisier three-point curvatures of very fine STL strips.
    for width in (0.002, 0.07):
        counts = tree.query_ball_point(keys, width, return_length=True)
        used = np.zeros(len(radius), dtype=bool)
        for index in np.argsort(-counts, kind='stable'):
            if used[index] or counts[index] < 2:
                continue
            near = tree.query_ball_point(keys[index], width)
            used[near] = True
            dirs = d[near]
            axis = (dirs * (height[near]**2 * np.sign(dirs @ d[index]))[:, None]).sum(axis=0)
            axis /= np.linalg.norm(axis)
            yield axis, np.median(center[near], axis=0), np.median(radius[near])


def _patch_cylinder(vertices, faces, normals, area, ids):
    """Fit a smooth patch without assuming any edge is an axial generator.

    The smallest normal covariance eigenvector initializes the axis. Fitting
    the vertices then removes the axis bias due to skew tessellation.
    """
    n = normals[ids]
    weights = area[ids]
    total = weights.sum()
    if total < 1e-12:
        return None
    values, vect = np.linalg.eigh(n.T @ (weights[:, None] * n) / total)
    if values[1] < 0.002 or values[0] > 0.18 * values[1]:
        return None
    axis = vect[:, 0]
    frame = _cyl_frame(axis)
    points = vertices[np.unique(faces[ids])]
    origin = points.mean(axis=0)
    delta = points - origin
    xy = delta @ frame.T
    fit = np.linalg.lstsq(np.column_stack((2*xy, np.ones(len(xy)))),
                          _cyl_dot(xy, xy), rcond=None)[0]
    r = np.sqrt(max(0., fit[2] + fit[:2] @ fit[:2]))
    if not (1.5 < r < 22.):
        return None
    center = fit[:2] @ frame
    radial = vertices[faces[ids,0]] - (center + origin)
    radial -= (radial @ axis)[:, None] * axis
    if np.sum(weights[_cyl_dot(radial, n) < -0.4*r]) < 0.8 * total:
        return None
    for _ in range(7):
        frame = _cyl_frame(axis)
        rel = delta - center
        xy = rel @ frame.T
        lengths = np.linalg.norm(xy, axis=1)
        err = lengths - r
        if np.max(np.abs(err)) < 0.00002:
            break
        radial = xy / np.maximum(lengths, 1e-12)[:,None]
        jac = np.column_stack((radial, (rel @ axis)[:, None]*radial, np.ones(len(xy))))
        change = np.linalg.lstsq(jac, err, rcond=1e-9)[0]
        if np.linalg.norm(change[2:4]) > 0.18 or abs(change[4]) > 2.:
            return None
        center += change[:2] @ frame
        axis += change[2:4] @ frame
        axis /= np.linalg.norm(axis)
        center -= (center @ axis) * axis
        r += change[4]
        if np.linalg.norm(change) < 1e-9:
            break
    if (not 1.7 < r < 20.1) or np.quantile(np.abs(err), 0.85) > 0.004:
        return None
    center += origin
    center -= (center @ axis) * axis
    return axis, center, r


def _smooth_seeds(vertices, faces, normals, area, adjacency):
    _, first, second = adjacency
    i, j = first // 3, second // 3
    smooth = _cyl_dot(normals[i], normals[j]) > np.cos(np.deg2rad(19.))
    i, j = i[smooth], j[smooth]
    graph = coo_matrix((np.ones(2*len(i), dtype=bool), (np.r_[i,j], np.r_[j,i])),
                        shape=(len(faces), len(faces))).tocsr()
    _, label = connected_components(graph, directed=False)
    order = np.argsort(label)
    marks = np.r_[0, 1 + np.flatnonzero(np.diff(label[order])), len(faces)]
    used = np.zeros(len(faces), dtype=bool)
    for lo, hi in zip(marks[:-1], marks[1:]):
        if hi - lo < 5:
            continue
        ids = order[lo:hi]
        candidate = _patch_cylinder(vertices, faces, normals, area, ids)
        if candidate is not None:
            yield candidate
            continue
        if hi - lo < 30:
            continue
        # Smooth fillets sometimes connect cylinders to other faces. Try
        # overlapping local neighborhoods within such a connected component.
        for i in ids:
            if used[i]:
                continue
            patch = {i}
            frontier = {i}
            for depth in range(9):
                nbrs = {v for u in frontier
                        for v in graph.indices[graph.indptr[u]:graph.indptr[u+1]]}
                frontier = nbrs.difference(patch)
                if not frontier:
                    break
                patch.update(frontier)
                if depth == 3:
                    used[list(patch)] = True
            used[i] = True
            if len(patch) < 8 or len(patch) == len(ids):
                continue
            candidate = _patch_cylinder(vertices, faces, normals, area,
                                        np.fromiter(patch, dtype=int))
            if candidate is not None:
                yield candidate


def _land_mask(vertices, faces, normals, axis, center, radius, tolerance, trial=False):
    delta = vertices - center
    radial = delta - (delta @ axis)[:, None] * axis
    lengths = np.linalg.norm(radial, axis=1)
    on = np.all(np.abs(lengths - radius)[faces] <= tolerance, axis=1)
    ids = np.flatnonzero(on)
    if len(ids):
        # Axial end faces can have all three vertices on a rim. Their normals
        # distinguish them from inward-facing wall.
        tri = faces[ids]
        inward = _cyl_dot(normals[ids], radial[tri[:, 0]]) < -0.05 * radius
        projected = radial[tri]
        other = np.roll(projected, 1, axis=1)
        dirs = projected / np.maximum(lengths[tri, None], 1e-12)
        cos_span = _cyl_dot(dirs, np.roll(dirs, 1, axis=1)).min(axis=1)
        # The angular extent is the largest pairwise separation for a small
        # facet. Its largest radial deviation is at the closest projected
        # point on its three edges (an axis-containing triangle fails span).
        side = other - projected
        t = np.clip(-_cyl_dot(projected, side) /
                    np.maximum(_cyl_dot(side, side), 1e-30), 0., 1.)
        closest = projected + t[:, :, None] * side
        sag = radius - np.sqrt(_cyl_dot(closest, closest).min(axis=1))
        max_angle = 17. if trial else 16.003
        max_sag = 0.05 + (tolerance if trial else 0.00008)
        ids = ids[inward & (cos_span >= np.cos(np.deg2rad(max_angle))) & (sag <= max_sag)]
    return ids


def _wall_faces(vertices, faces, normals, axis, center, radius, seeds):
    """Collect trimmed triangles adjacent to faces on the nominal cylinder.

    Curve-intersection vertices from the CAD tessellator can sit inward of
    the exact radius. Grow across the wall, with the facet chord and dihedral
    constraints preventing propagation through a shoulder or an end cap.
    """
    near = _land_mask(vertices, faces, normals, axis, center, radius, 0.0501)
    if len(seeds) == len(near):
        return near
    # Include detached patches with tight vertex evidence in seeds before the
    # flood so that connectivity is never a condition on the entire land.
    on = np.zeros(len(faces), bool)
    on[seeds] = True
    near_faces = faces[near]
    edges = np.stack((near_faces, np.roll(near_faces, -1, axis=1)), axis=-1).reshape(-1,2)
    edges.sort(axis=1)
    keys = edges[:, 0]*len(vertices) + edges[:, 1]
    order = np.argsort(keys)
    equal = keys[order[:-1]] == keys[order[1:]]
    left = near[order[:-1][equal] // 3]
    right = near[order[1:][equal] // 3]
    adjacent = _cyl_dot(normals[left], normals[right])
    # Interior and singly trimmed wall facets have small vertex errors.
    # A boolean cut through an existing STL chord can leave larger errors,
    # up to the full sag. In that case extend the same chord plane only.
    # Traversing arbitrary small dihedrals at the full sag tolerance could
    # absorb nearby spline surfaces, e.g. the curved outline of lettering.
    near_nominal = np.zeros(len(faces), bool)
    near_nominal[_land_mask(vertices, faces, normals, axis, center, radius, 0.006)] = True
    coplanar = adjacent > np.cos(np.deg2rad(0.06))
    gentle = adjacent > np.cos(np.deg2rad(20.))
    step_left = gentle & (near_nominal[left] | coplanar)
    step_right = gentle & (near_nominal[right] | coplanar)
    while True:
        reach_left = step_left & on[right] & ~on[left]
        reach_right = step_right & on[left] & ~on[right]
        if not np.any(reach_left) and not np.any(reach_right):
            break
        on[left[reach_left]] = True
        on[right[reach_right]] = True
    return np.flatnonzero(on)


def _refit_land(vertices, faces, ids, axis, center, radius):
    """Robust orthogonal-distance fit, using distinct vertices in a local frame.

    First refine the circle perpendicular to the generator or rim normal. Then
    allow small axis corrections. Limiting the influence of trimmed vertices
    keeps the unconstrained fit from tilting towards chord-intersection points.
    """
    points = vertices[np.unique(faces[ids])]
    origin = points.mean(axis=0)
    delta = points - origin
    axis = axis.copy()
    center = center - origin
    center -= (center @ axis) * axis
    for fixed in (True, False):
        for _ in range(12):
            frame = _cyl_frame(axis)
            rel = delta - center
            xy = rel @ frame.T
            lengths = np.linalg.norm(xy, axis=1)
            residual = lengths - radius
            normals = xy / np.maximum(lengths, 1e-12)[:, None]
            scale = max(2.5e-5, 0.85 * np.median(np.abs(residual)))
            weights = 1. / np.sqrt(1. + (residual/scale)**2)
            columns = [normals, np.ones(len(points))]
            if not fixed:
                columns.append((rel @ axis)[:, None] * normals)
            jac = np.column_stack(columns) * weights[:, None]
            change = np.linalg.lstsq(jac, residual * weights, rcond=1e-10)[0]
            if not fixed and np.linalg.norm(change[3:]) > 0.015:
                break
            center += change[:2] @ frame
            radius += change[2]
            if not fixed:
                axis += change[3:] @ frame
                axis /= np.linalg.norm(axis)
                center -= (center @ axis) * axis
            if np.linalg.norm(change) < 3e-10:
                break
    center += origin
    center -= (center @ axis) * axis
    return axis, center, radius


def _occupied_angle(vertices, faces, ids, axis, center):
    frame = _cyl_frame(axis)
    xy = (vertices - center) @ frame.T
    a = np.arctan2(xy[:, 1], xy[:, 0])[faces[ids]]
    local = (a - a[:, :1] + np.pi) % (2*np.pi) - np.pi
    lo = (a[:, 0] + np.min(local, axis=1)) % (2*np.pi)
    hi = lo + np.ptp(local, axis=1)
    spans = sorted(zip(np.r_[lo, np.zeros(np.count_nonzero(hi > 2*np.pi))],
                       np.r_[np.minimum(hi, 2*np.pi), hi[hi > 2*np.pi] - 2*np.pi]))
    length = 0.
    end = 0.
    for a, b in spans:
        length += max(0., b - max(a, end))
        end = max(b, end)
    return length


def find_lands(triangles):
    """Detect inward, nominal cylindrical lands in STL coordinates."""
    vertices, faces = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    faces = faces.reshape(-1, 3)
    origin = vertices.mean(axis=0)
    vertices -= origin
    normals = np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                       vertices[faces[:, 2]] - vertices[faces[:, 0]])
    area = np.linalg.norm(normals, axis=1)
    normals /= np.maximum(1e-30, area)[:, None]
    adjacency = _land_adjacency(faces, len(vertices))
    edges = _crease_circles(vertices, faces, normals, adjacency)
    seeds = list(_crease_seeds(*edges))
    seeds.extend(rim_seeds(vertices, normals, adjacency))
    seeds.extend(_smooth_seeds(vertices, faces, normals, area, adjacency))
    lands = []
    for d, c, r in seeds:
        if any(np.linalg.norm(np.cross(d, old_d)) < 4e-6 and
               np.linalg.norm(c - old_c) < 3e-4 and abs(r - old_r) < 3e-4
               for old_c, old_d, old_r, _ in lands):
            continue
        # Use the most accurate seed evidence first so two close radii do not
        # get averaged into a single unsupported cylinder.
        for width in (0.00022, 0.003, 0.035):
            ids = _land_mask(vertices, faces, normals, d, c, r, width, trial=True)
            if len(ids) >= 5:
                break
        if len(ids) < 5:
            continue
        d, c, r = _refit_land(vertices, faces, ids, d, c, r)
        ids = _land_mask(vertices, faces, normals, d, c, r, 0.00022)
        if len(ids) < 5:
            ids = _land_mask(vertices, faces, normals, d, c, r, 0.001)
        if len(ids) < 5:
            continue
        d, c, r = _refit_land(vertices, faces, ids, d, c, r)
        if not 1.75 - 3e-5 <= r <= 20 + 3e-5:
            continue
        ids = _land_mask(vertices, faces, normals, d, c, r, 0.00022)
        if len(ids) < 5:
            continue
        ids = _wall_faces(vertices, faces, normals, d, c, r, ids)
        if _occupied_angle(vertices, faces, ids, d, c) < np.deg2rad(160. - 0.002):
            continue
        if any(np.array_equal(ids, old_ids) or
               (np.linalg.norm(np.cross(d, old_d)) < 4e-6 and
                np.linalg.norm(c - old_c) < 3e-4 and abs(r - old_r) < 3e-4)
               for old_c, old_d, old_r, old_ids in lands):
            continue
        d, c, r = _refit_land(vertices, faces, ids, d, c, r)
        lands.append((c, d, float(r), ids))
    return [(c + origin - ((c + origin) @ d) * d, d, r, ids)
            for c, d, r, ids in lands]


def _rim_triplets(vertices, normals, adjacency):
    """Find consecutive edges that turn gradually along mesh creases."""
    sides, first, second = adjacency
    # Small but distinct dihedrals retain the rims at shallow chamfers.
    crease = _cyl_dot(normals[first // 3], normals[second // 3]) < np.cos(np.deg2rad(0.4))
    edges = sides[first[crease]]
    if len(edges) == 0:
        return np.empty((0, 3), dtype=int)
    directed = np.concatenate((edges, edges[:, ::-1]))
    directed = directed[np.argsort(directed[:, 0])]
    counts = np.bincount(directed[:, 0], minlength=len(vertices))
    starts = counts.cumsum() - counts
    triples = []
    # Group vertices by valence to compare their incident directions in batches.
    for degree in np.unique(counts[counts >= 2]):
        ids = np.flatnonzero(counts == degree)
        i, j = np.triu_indices(degree, 1)
        neighbors = directed[starts[ids, None] + np.arange(degree), 1]
        directions = vertices[neighbors] - vertices[ids, None]
        directions /= np.maximum(np.linalg.norm(directions, axis=2, keepdims=True), 1e-30)
        dot = _cyl_dot(directions[:, i], directions[:, j])
        row, pair = np.where((dot < -np.cos(np.deg2rad(16.1))) &
                             (dot > -np.cos(np.deg2rad(0.04))))
        triples.append(np.column_stack((neighbors[row, i[pair]], ids[row], neighbors[row, j[pair]])))
    return np.concatenate(triples) if triples else np.empty((0, 3), dtype=int)


def _rim_circles(vertices, triples):
    p, q, s = vertices[triples].transpose(1, 0, 2)
    x, y = p - q, s - q
    normal = np.cross(x, y)
    length = np.linalg.norm(normal, axis=1)
    usable = length > 1e-12
    normal, length, x, y, q = (a[usable] for a in (normal, length, x, y, q))
    # Circumcenter in the plane of each triple, based at its middle vertex.
    offset = (np.cross(normal, x) * _cyl_dot(y, y)[:, None] +
              np.cross(y, normal) * _cyl_dot(x, x)[:, None]) / (2 * length[:, None]**2)
    radius = np.linalg.norm(offset, axis=1)
    keep = (radius > 1.70) & (radius < 20.10)
    return (normal / length[:, None])[keep], (q + offset)[keep], radius[keep]


def rim_seeds(vertices, normals, adjacency):
    """Estimate cylinders from circular rim chains, including staggered rims."""
    short = _rim_triplets(vertices, normals, adjacency)
    if len(short) == 0:
        return []
    # Wide triples suppress the sensitivity of three-point curvature to STL
    # rounding. Retain the original triples as well for interrupted rims.
    following = {}
    for a, b, c in short:
        for edge, end in (((a, b), c), ((c, b), a)):
            following[edge] = end if edge not in following else None
    wide = short.copy()
    for i, (a, b, c) in enumerate(short):
        for k, first in ((0, a), (2, c)):
            v, w = b, first
            for _ in range(3):
                nxt = following.get((v, w))
                if nxt is None:
                    break
                v, w = w, nxt
            wide[i, k] = w
    axes, centers, radii = _rim_circles(vertices, np.concatenate((short, wide)))
    if len(axes) == 0:
        return []
    centers -= _cyl_dot(centers, axes)[:, None] * axes
    ii, jj = np.triu_indices(3)
    features = np.column_stack((15 * axes[:, ii] * axes[:, jj], centers, radii))
    tree = cKDTree(features)
    counts = tree.query_ball_point(features, 0.018, return_length=True)
    occupied = np.zeros(len(axes), dtype=bool)
    seeds = []
    for i in np.argsort(-counts):
        if occupied[i] or counts[i] < 2:
            continue
        ids = tree.query_ball_point(features[i], 0.013)
        occupied[ids] = True
        directions = axes[ids] * np.sign(axes[ids] @ axes[i])[:, None]
        axis = directions.mean(axis=0)
        axis /= np.linalg.norm(axis)
        seeds.append((axis, np.median(centers[ids], axis=0), np.median(radii[ids])))
    return seeds


@njit
def _disk_triangle_extrema(tri, cx, cy, radius):
    """Extrema of z on the part of one triangle inside the vertical cylinder."""
    low = np.inf
    high = -np.inf
    tol = 2e-12  # roundoff in distances, mm
    # The endpoints of each clipped edge include all in-disk vertices.
    for i in range(3):
        j = (i + 1) % 3
        x = tri[i, 0] - cx
        y = tri[i, 1] - cy
        dx = tri[j, 0] - tri[i, 0]
        dy = tri[j, 1] - tri[i, 1]
        dz = tri[j, 2] - tri[i, 2]
        length2 = dx * dx + dy * dy
        if length2 <= tol * tol:
            if x * x + y * y <= (radius + tol) ** 2:
                low = min(low, tri[i, 2], tri[j, 2])
                high = max(high, tri[i, 2], tri[j, 2])
            continue
        length = math.sqrt(length2)
        distance = abs(x * dy - y * dx) / length
        if distance > radius + tol:
            continue
        middle = -(x * dx + y * dy) / length2
        half = math.sqrt(max(0., (radius - distance) * (radius + distance))) / length
        a = max(0., middle - half)
        b = min(1., middle + half)
        if a <= b + tol / length:
            za = tri[i, 2] + dz * min(a, 1.)
            zb = tri[i, 2] + dz * max(b, 0.)
            low = min(low, za, zb)
            high = max(high, za, zb)

    # On a nonvertical facet z is affine in (x,y). On an arc of the
    # circle its only stationary points are along the projected gradient.
    ax = tri[1, 0] - tri[0, 0]
    ay = tri[1, 1] - tri[0, 1]
    az = tri[1, 2] - tri[0, 2]
    bx = tri[2, 0] - tri[0, 0]
    by = tri[2, 1] - tri[0, 1]
    bz = tri[2, 2] - tri[0, 2]
    det = ax * by - ay * bx
    if det == 0.:
        return low, high
    gx = az * by - ay * bz
    gy = ax * bz - az * bx
    scale = math.hypot(gx, gy)
    if scale != 0.:
        gx = radius * gx / scale
        gy = radius * gy / scale
    for sign in (-1., 1.):
        x = (cx - tri[0, 0]) + sign * gx
        y = (cy - tri[0, 1]) + sign * gy
        u = (x * by - y * bx) / det
        v = (ax * y - ay * x) / det
        if u >= -1e-13 and v >= -1e-13 and u + v <= 1. + 1e-13:
            z = tri[0, 2] + u * az + v * bz
            low = min(low, z)
            high = max(high, z)
    return low, high


@njit
def sweep_extrema(local_triangles, centers):
    """Return the first surface height seen from -axis and +axis, respectively."""
    radius = 1.5
    margin = radius + 2e-12
    bounds = np.empty((len(local_triangles), 4))
    for j, tri in enumerate(local_triangles):
        bounds[j, 0] = min(tri[0, 0], tri[1, 0], tri[2, 0]) - margin
        bounds[j, 1] = max(tri[0, 0], tri[1, 0], tri[2, 0]) + margin
        bounds[j, 2] = min(tri[0, 1], tri[1, 1], tri[2, 1]) - margin
        bounds[j, 3] = max(tri[0, 1], tri[1, 1], tri[2, 1]) + margin
    lows = np.full(len(centers), np.inf)
    highs = np.full(len(centers), -np.inf)
    for i, (cx, cy) in enumerate(centers):
        for j in range(len(local_triangles)):
            if cx < bounds[j, 0] or cx > bounds[j, 1]:
                continue
            if cy < bounds[j, 2] or cy > bounds[j, 3]:
                continue
            low, high = _disk_triangle_extrema(local_triangles[j], cx, cy, radius)
            lows[i] = min(lows[i], low)
            highs[i] = max(highs[i], high)
    return lows, highs


def transverse_clock(triangles, c, d):
    """Choose an axis sense and clock using the geometry of the component."""
    vertices = triangles.reshape(-1, 3)
    midpoint = vertices.mean(axis=0)
    relative_station = (vertices - midpoint) @ d
    if np.mean(relative_station**3) < 0:
        d = -d
    radial = midpoint - c
    radial -= d * np.dot(radial, d)
    if np.linalg.norm(radial) < 1e-4:
        displacement = vertices - c
        vectors = displacement - (displacement @ d)[:, None] * d
        radial = vectors[np.argmax(np.linalg.norm(vectors, axis=1))]
    radial /= np.linalg.norm(radial)
    # Offset from the part symmetry lines to avoid sampling exactly at trims.
    angle = np.deg2rad(0.371)
    e = np.cos(angle) * radial + np.sin(angle) * np.cross(d, radial)
    return d, e


def wall_section(wall_triangles, n):
    """Intersect the wall facets with the radial half-plane at this clock."""
    xy = wall_triangles[:, :, :2]
    side = xy @ np.array([-n[1], n[0]])
    side[np.abs(side) < 2e-11] = 0.0
    crosses = ((side.min(axis=1) <= 0) & (side.max(axis=1) >= 0) &
               ((xy.sum(axis=1) @ n) > 0))
    side = side[crosses]
    if not len(side):
        return []
    height = wall_triangles[crosses, :, 2]
    low = np.min(np.where(side == 0, height, np.inf), axis=1)
    high = np.max(np.where(side == 0, height, -np.inf), axis=1)
    for i, j in ((0, 1), (1, 2), (2, 0)):
        edge = (side[:, i] * side[:, j]) < 0
        t = side[edge, i] / (side[edge, i] - side[edge, j])
        z = height[edge, i] + t * (height[edge, j] - height[edge, i])
        low[edge] = np.minimum(low[edge], z)
        high[edge] = np.maximum(high[edge], z)

    intervals = []
    for i in np.argsort(low):
        a, b = float(low[i]), float(high[i])
        if b <= a:
            continue
        if intervals and a <= intervals[-1][1] + 2e-8:
            intervals[-1][1] = max(intervals[-1][1], b)
        else:
            intervals.append([a, b])
    return intervals


def survey_land(triangles, c, d, r, face_ids):
    d, e = transverse_clock(triangles, c, d)
    basis = np.column_stack((e, np.cross(d, e), d))
    # c is perpendicular to d, so local z is the absolute axial station h.
    local = np.ascontiguousarray((triangles - c) @ basis)
    wall = local[face_ids]
    theta = np.deg2rad(2.0 * np.arange(180))
    directions = np.column_stack((np.cos(theta), np.sin(theta)))
    lows, highs = sweep_extrema(local, (r - 1.70) * directions)
    rows = []
    for radial, low, high in zip(directions, lows, highs):
        section = wall_section(wall, radial)
        # From either infinity the first point of a closed, bounded solid lies
        # on its boundary. The window stops two millimetres behind that point.
        plus_limit = high + 2.0
        minus_limit = low - 2.0
        plus = [[max(a, plus_limit), b] for a, b in section if b > plus_limit]
        minus = [[a, min(b, minus_limit)] for a, b in section if a < minus_limit]
        rows.append([section, plus, minus])
    return dict(c=c.tolist(), d=d.tolist(), e=e.tolist(), r=float(r), rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    try:
        triangles = read_triangles(args.mesh)
        bores = [survey_land(triangles, *land) for land in find_lands(triangles)]
        args.out.mkdir(parents=True, exist_ok=True)
        fd, path = tempfile.mkstemp(dir=args.out, prefix=".survey-")
        tmp = Path(path)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as stream:
                json.dump({"bores": bores}, stream, allow_nan=False, separators=(",", ":"))
                stream.write("\n")
            os.replace(tmp, args.out / "survey.json")
        finally:
            tmp.unlink(missing_ok=True)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{parser.prog}: {exc}\n")


if __name__ == "__main__":
    main()
BORE_ACCESS_PY
chmod 755 "$out/bore_access.py"

/opt/conda/bin/python "$out/bore_access.py" --mesh "/app/inputs/Scope Base v33.stl" --out "$out/base"
/opt/conda/bin/python "$out/bore_access.py" --mesh "/app/inputs/Scope Screen Leg Mount v33.stl" --out "$out/screen-leg"
/opt/conda/bin/python "$out/bore_access.py" --mesh "/app/inputs/Scope Top Stage v33.stl" --out "$out/top-stage"
