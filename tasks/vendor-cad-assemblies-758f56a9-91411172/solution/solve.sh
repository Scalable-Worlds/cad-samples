#!/bin/bash
set -euo pipefail
mkdir -p /app/out
rm -f /app/out/pull.py /app/out/axial.ply /app/out/side.ply /app/out/inclined.ply
cat > /app/out/pull.py <<'PULL_PYTHON'
#!/opt/conda/bin/python
"""Construct the polyhedral steel trapped by a two-sided mould pull.

All incidence, depth ordering and interpolation use rational coordinates.
Spatial trees in millimetres accelerate bounding-box searches. They never
decide incidence. The final PLY stores file-frame XYZ in binary64.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import defaultdict, deque
from functools import cmp_to_key
from itertools import islice
from math import gcd
from pathlib import Path
from tempfile import TemporaryDirectory

import gmpy2
import numpy as np
import shapely
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

Q = gmpy2.mpq


def stock(path):
    mesh = trimesh.load_mesh(path, process=False)
    _, first, ids = np.unique(np.round(mesh.vertices, 5), axis=0,
                              return_index=True, return_inverse=True)
    # Use the first stored coordinate for each weld key, preserving the STL
    # facet positions. A common binary denominator makes predicates integral.
    raw = mesh.vertices[first]
    ratios = [[float(t).as_integer_ratio() for t in p] for p in raw]
    scale = max(b for row in ratios for a, b in row)
    v = [tuple(a * (scale // b) for a, b in row) for row in ratios]
    faces = []
    seen = set()
    for f in ids[mesh.faces].tolist():
        key = tuple(sorted(f))
        if key in seen:
            continue
        seen.add(key)
        a, b, c = (v[i] for i in f)
        s = [y-x for x, y in zip(a, b)]
        t = [y-x for x, y in zip(a, c)]
        if any(s[i]*t[j] != s[j]*t[i] for i, j in ((0, 1), (1, 2), (2, 0))):
            faces.append(f)
    return v, faces, scale


def determinant(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


def pull_frame(vertices, direction, scale):
    """Use two integral projection coordinates and depth along a file axis."""
    k = max(range(3), key=lambda axis: abs(direction[axis]))
    i, j = (k + 1) % 3, (k + 2) % 3
    w = direction[k]
    origin = vertices[0]
    points = []
    for vertex in vertices:
        t = tuple(a - b for a, b in zip(vertex, origin))
        points.append((w*t[i] - direction[i]*t[k], w*t[j] - direction[j]*t[k], t[k]))

    def to_world(q):
        x, y, h = q
        world = [0.0] * 3
        world[k] = float(Q(h + origin[k], scale))
        world[i] = float(Q(x + direction[i]*h + w*origin[i], w*scale))
        world[j] = float(Q(y + direction[j]*h + w*origin[j], w*scale))
        return world

    return points, to_world


def plane(tri):
    p, q, r = tri
    dx, dy, dz = [b-a for a, b in zip(p, q)]
    ex, ey, ez = [b-a for a, b in zip(p, r)]
    z = dx*ey-ex*dy
    a = dz*ey-ez*dy
    b = dx*ez-ex*dz
    c = z*p[2]-a*p[0]-b*p[1]
    return a, b, c, z


def elevation(s, x, y):
    a, b, c, z = s
    return Q(a*x+b*y+c, z)


def merged_lines(triangles):
    # A supporting line is y=(a*x+c)/b with b>0. First combine overlapping
    # collinear segments. Separate runs on the same line can enter and leave
    # the sweep independently.
    intervals = defaultdict(list)
    vertical = defaultdict(list)
    for t in triangles:
        for p, q in zip(t, (t[1], t[2], t[0])):
            if p[0] == q[0]:
                vertical[p[0]].append(sorted((p[1], q[1])))
                continue
            if p[0] > q[0]:
                p, q = q, p
            b, a = q[0]-p[0], q[1]-p[1]
            common = gcd(b, a)
            a, b = a // common, b // common
            c = b*p[1]-a*p[0]
            intervals[a, b, c].append((p[0], q[0]))
    lines = list(intervals)
    runs = []
    for label, line in enumerate(lines):
        spans = sorted(intervals[line])
        l, r = spans[0]
        for p, q in spans[1:]:
            if p <= r:
                r = max(r, q)
            else:
                runs.append((label, l, r))
                l, r = p, q
        runs.append((label, l, r))
    return lines, runs, vertical


def sweep_events(lines, runs, scale):
    # Each event holds lines to remove and then insert in the vertical order.
    events = defaultdict(lambda: (set(), set()))
    boxes = []
    for i, l, r in runs:
        events[l][1].add(i)
        events[r][0].add(i)
        a, b, c = lines[i]
        y0, y1 = (a*l+c) // b, (a*r+c) // b
        boxes.append((l / scale, min(y0, y1) / scale, r / scale, max(y0, y1) / scale))
    boxes = shapely.box(*np.asarray(boxes, dtype=float).T)
    tree = shapely.STRtree(boxes)
    for start in range(0, len(runs), 128):
        js, ks = tree.query(boxes[start:start+128])
        for j, k in zip((js+start).tolist(), ks.tolist()):
            if j >= k:
                continue
            i, l, r = runs[j]
            m, s, t = runs[k]
            a, b, c = lines[i]
            e, f, g = lines[m]
            den = a*f-e*b
            if not den:
                continue
            x = Q(g*b-c*f, den)
            if max(l, s) <= x <= min(r, t):
                rem, add = events[x]
                if l < x < r:
                    rem.add(i)
                    add.add(i)
                if s < x < t:
                    rem.add(m)
                    add.add(m)
    return events


def trapezoids(lines, events, vertical):
    """Vertical decomposition of the projected straight-line arrangement."""
    active = []
    live = {}
    for x in sorted(events.keys() | vertical.keys()):
        remove, insert = events[x]
        def close(p, q):
            if p is not None and q is not None:
                left = live.pop((p, q))
                if left < x:
                    return [(left, x, p, q)]
            return []

        def begin(p, q):
            if p is not None and q is not None:
                live[p, q] = x

        # Stock edges at constant abscissa can separate different stacks even
        # when both delimiting nonvertical lines continue through this event.
        for lower, upper in vertical.get(x, ()):
            lo, hi = 0, len(active)
            while lo < hi:
                mid = (lo+hi) // 2
                if ordinate(lines, active[mid], x) <= lower:
                    lo = mid+1
                else:
                    hi = mid
            for pos in range(max(lo, 1), len(active)):
                p, q = active[pos-1], active[pos]
                if ordinate(lines, p, x) >= upper:
                    break
                yield from close(p, q)
                begin(p, q)
        for i in sorted(remove):
            pos = active.index(i)
            low = active[pos-1] if pos else None
            high = active[pos+1] if pos+1 < len(active) else None
            yield from close(low, i)
            yield from close(i, high)
            active.pop(pos)
            begin(low, high)
        for i in sorted(insert):
            a, b, c = lines[i]
            lo, hi = 0, len(active)
            while lo < hi:
                mid = (lo+hi) // 2
                e, f, g = lines[active[mid]]
                slope = e*b-a*f
                sign = slope*x+g*b-c*f
                # At a tied ordinate the open slab to the right determines order.
                if sign < 0 or (sign == 0 and slope < 0):
                    lo = mid+1
                else:
                    hi = mid
            low = active[lo-1] if lo else None
            high = active[lo] if lo < len(active) else None
            yield from close(low, high)
            active.insert(lo, i)
            begin(low, i)
            begin(i, high)


def ordinate(lines, i, x):
    a, b, c = lines[i]
    return Q(a*x+c, b)


def occupied_cells(triangles, planes, lines, events, vertical, scale):
    """Classify one strict interior ray per trapezoid; planes cannot swap there."""
    t = np.array([[(p[0] / scale, p[1] / scale) for p in tri] for tri in triangles], dtype=float)
    lo, hi = t.min(axis=1), t.max(axis=1)
    tree = shapely.STRtree(shapely.box(lo[:, 0], lo[:, 1], hi[:, 0], hi[:, 1]))
    slabs = trapezoids(lines, events, vertical)
    while batch := list(islice(slabs, 512)):
        rays = []
        for l, r, i, j in batch:
            x = Q(l+r, 2)
            y = (ordinate(lines, i, x)+ordinate(lines, j, x))/2
            rays.append((x, y))
        arr = shapely.points([(float(x/scale), float(y/scale)) for x, y in rays])
        qi, fi = tree.query(arr)
        hits = [[] for _ in batch]
        for ray, top in zip(qi.tolist(), fi.tolist()):
            p = rays[ray]
            tri = triangles[top]
            # Cross multiplication gives the strict interior predicate even for
            # a microscopic cell whose midpoint rounds onto an edge in double.
            if all(determinant(a, b, p) > 0 for a, b in zip(tri, (tri[1], tri[2], tri[0]))):
                hits[ray].append((elevation(planes[top], *p), top))
        for cell, h in zip(batch, hits):
            if len(h) % 2:
                raise ValueError("Unpaired stock crossing in an open projection cell")
            if len(h) > 2:
                h.sort()
                spans = [(a[1], b[1]) for a, b in zip(h[1::2], h[2::2]) if b[0] > a[0]]
                if spans:
                    yield (*cell, spans)


def air_boundary(p, faces, scale, to_world):
    triangles = []
    # A facet parallel to the pull has no open footprint. Its exposed steel
    # walls arise from the difference of the spans on the two sides.
    for f in faces:
        a, b, c = [p[i] for i in f]
        orient = determinant(a, b, c)
        if orient:
            triangles.append((a, b, c) if orient > 0 else (a, c, b))
    planes = [plane(t) for t in triangles]
    lines, runs, vertical = merged_lines(triangles)
    events = sweep_events(lines, runs, scale)
    cells = list(occupied_cells(triangles, planes, lines, events, vertical, scale))

    # Node the boundaries of the occupied trapezoids, including the vertical
    # edges introduced by the sweep. Both neighbors must use the same sequence
    # of vertices along a shared edge, even if one is much longer.
    on_line, on_vertical = defaultdict(set), defaultdict(set)
    for l, r, a, b, spans in cells:
        for line in (a, b):
            on_line[line].update((l, r))
            for x in (l, r):
                on_vertical[x].add(ordinate(lines, line, x))
    for index in (on_line, on_vertical):
        for key in index:
            index[key] = sorted(index[key])
    xy = []
    index = {}

    def node(x, y):
        if (x, y) not in index:
            index[x, y] = len(xy)
            xy.append((x, y))
        return index[x, y]

    def pieces():
        for l, r, a, b, spans in cells:
            la, ra = ordinate(lines, a, l), ordinate(lines, a, r)
            lb, rb = ordinate(lines, b, l), ordinate(lines, b, r)
            lower, upper = on_line[a], on_line[b]
            left, right = on_vertical[l], on_vertical[r]
            ring = [node(x, ordinate(lines, a, x)) for x in lower[bisect_left(lower, l):bisect_right(lower, r)]]
            ring.extend(node(r, y) for y in right[bisect_right(right, ra):bisect_left(right, rb)])
            ring.extend(node(x, ordinate(lines, b, x)) for x in reversed(upper[bisect_left(upper, l):bisect_right(upper, r)]))
            ring.extend(node(l, y) for y in reversed(left[bisect_right(left, la):bisect_left(left, lb)]))
            if ring[0] == ring[-1]:
                ring.pop()
            if ra == rb:
                k = ring.index(node(r, ra))
                ring.pop(k)
            if len(ring) == 3:
                facets = [tuple(ring)]
            else:
                m = node(Q(l+r, 2), (la+lb+ra+rb)/4)
                facets = [(i, j, m) for i, j in zip(ring, ring[1:]+ring[:1])]
            for low, high in spans:
                ends = {i:(elevation(planes[low], *xy[i]), elevation(planes[high], *xy[i]))
                        for i in dict.fromkeys(i for facet in facets for i in facet)}
                for i, j, k in facets:
                    lo, hi = zip(ends[i], ends[j], ends[k])
                    yield i, j, k, lo, hi

    return stitch(xy, pieces(), to_world)


def _exact_planes(coords, faces):
    """Label each supporting plane and its signed invertible 2D projection."""
    table = {}
    axes = []
    plane_id = []
    sign = []
    for tri in faces:
        p, q, r = (coords[k] for k in tri)
        ux, uy, uz = (b-a for a, b in zip(p, q))
        vx, vy, vz = (b-a for a, b in zip(p, r))
        n = (uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx)
        axis = 2 if n[2] else 0 if n[0] else 1
        a, b = (axis+1) % 3, (axis+2) % 3
        l = gmpy2.mpq(n[a], n[axis])
        m = gmpy2.mpq(n[b], n[axis])
        key = (axis, l, m, p[axis]+l*p[a]+m*p[b])
        label = table.get(key)
        if label is None:
            label = table[key] = len(axes)
            axes.append((a, b))
        plane_id.append(label)
        sign.append(1 if n[axis] > 0 else -1)
    return axes, plane_id, sign


def simplify_faces(coords, facets):
    """Exact link-preserving contractions on planar patches and straight creases."""
    faces = np.asarray(facets, dtype=np.int64).tolist()
    axes, planes, signs = _exact_planes(coords, faces)
    incident = [set() for _ in coords]
    for f, vertices in enumerate(faces):
        for i in vertices:
            incident[i].add(f)
    queue = deque(range(len(coords)))
    queued = bytearray([1]) * len(coords)
    while queue:
        src = queue.popleft()
        queued[src] = 0
        star = incident[src]
        # Leaving a high-valence corner untouched bounds the work of individual
        # candidates. This changes only how much redundant triangulation remains.
        if not 3 <= len(star) <= 64:
            continue
        ring = {}
        labels = set()
        for f in star:
            t = faces[f]
            j = t.index(src)
            a, b = t[(j+1) % 3], t[(j+2) % 3]
            if a in ring:
                break
            ring[a] = (b, f)
            labels.add(planes[f])
            if len(labels) > 2:
                break
        else:
            # Only regular, single-cycle face stars are eligible. This also
            # protects point and edge contacts before vertex-link separation.
            nxt = next(iter(ring))
            start = nxt
            order = []
            while nxt in ring and len(order) < len(ring):
                order.append(nxt)
                nxt = ring[nxt][0]
                if nxt == start:
                    break
            if nxt != start or len(order) != len(ring):
                continue
            degree = len(order)
            candidates = []
            for k, dst in enumerate(order):
                prev = order[k-1]
                side1, side2 = ring[prev][1], ring[dst][1]
                if len(labels) == 2 and planes[side1] == planes[side2]:
                    continue
                if len(incident[dst]) > 64:
                    continue
                candidates.append((len(incident[dst]), dst, prev, side1, side2))
            for _, dst, prev, side1, side2 in sorted(candidates):
                # An edge contraction must not identify neighboring edges or
                # faces, other than the two removed triangles and their sides.
                common = set()
                for f in incident[dst]:
                    common.update(v for v in faces[f] if v in ring)
                if common != {dst, prev, ring[dst][0]}:
                    continue
                # A four-vertex tetrahedron also has the entire opposing edge
                # in both links; contracting it creates coincident triangles.
                if degree == 3 and any(prev in faces[f] and ring[dst][0] in faces[f] for f in incident[dst]):
                    continue
                p = coords[dst]
                for a, (b, f) in ring.items():
                    if f in (side1, side2):
                        continue
                    x, y = axes[planes[f]]
                    q, r = coords[a], coords[b]
                    cross = (q[x]-p[x])*(r[y]-p[y])-(q[y]-p[y])*(r[x]-p[x])
                    if (cross if signs[f] > 0 else -cross) <= 0:
                        break
                else:
                    for f in star:
                        if f == side1 or f == side2:
                            for v in faces[f]:
                                if v != src:
                                    incident[v].remove(f)
                            faces[f] = None
                        else:
                            face = faces[f]
                            face[face.index(src)] = dst
                            incident[dst].add(f)
                    star.clear()
                    for near in order:
                        if not queued[near]:
                            queue.append(near)
                            queued[near] = 1
                    break
    return np.array([t for t in faces if t is not None], dtype=np.int64).reshape(-1, 3)


def _pair_junction(edge, halfedges, corner_vertices, coords):
    """Match the two faces that bound each local occupied angular sector."""
    p, q = (coords[k] for k in edge)
    axis = tuple(b-a for a, b in zip(p, q))
    ax, ay, az = axis
    ref = (-ay, ax, 0) if ax or ay else (1, 0, 0)
    rx, ry, rz = ref
    up = (ay*rz-az*ry, az*rx-ax*rz, ax*ry-ay*rx)
    rays = []
    for h in halfedges:
        base, j = divmod(int(h), 3)
        t = coords[corner_vertices[3*base+(j+2) % 3]]
        v = tuple(b-a for a, b in zip(p, t))
        a = sum(k*l for k, l in zip(v, ref))
        b = sum(k*l for k, l in zip(v, up))
        half = 0 if b > 0 or (b == 0 and a >= 0) else 1
        rays.append((half, a, b, int(h)))

    def compare(i, j):
        if i[0] != j[0]:
            return i[0]-j[0]
        orient = i[1]*j[2] - i[2]*j[1]
        if orient:
            return -1 if orient > 0 else 1
        return i[3]-j[3]

    rays.sort(key=cmp_to_key(compare))
    pairs = []
    for k, ray in enumerate(rays):
        h = ray[3]
        if corner_vertices[h] == edge[0]:
            other = rays[k-1][3]
            if corner_vertices[other] != edge[1]:
                raise ValueError("Inconsistent angular order along a boundary junction")
            pairs.append((h, other))
    if len(pairs)*2 != len(halfedges):
        raise ValueError("Open boundary at a junction")
    return pairs


def _normalize_boundary(coords, triangles, to_world):
    """Give touching sheets independent indices by following their face links."""
    facets = np.asarray(triangles, dtype=np.int64)
    # Two gaps can meet on a whole cap if the intervening stock has tapered
    # to zero thickness. Such coincident, opposite faces are interior.
    keys = np.sort(facets, axis=1)
    order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
    same = np.all(keys[order[:-1]] == keys[order[1:]], axis=1)
    if np.any(same):
        keep = np.ones(len(facets), dtype=bool)
        keep[order[:-1][same]] = keep[order[1:][same]] = False
        facets = facets[keep]
    if len(facets) == 0:
        return np.empty((0, 3)), np.empty((0, 3), dtype=np.int64)
    facets = simplify_faces(coords, facets)
    corners = facets.reshape(-1)
    next_corner = np.arange(len(corners), dtype=np.int64).reshape(-1, 3)[:, [1, 2, 0]].ravel()
    end = corners[next_corner]
    low, high = np.minimum(corners, end), np.maximum(corners, end)
    order = np.lexsort((high, low))
    cuts = np.r_[0, np.flatnonzero((np.diff(low[order]) != 0) |
                                   (np.diff(high[order]) != 0)) + 1, len(order)]
    counts = np.diff(cuts)
    ordinary = cuts[:-1][counts == 2]
    first = order[ordinary]
    second = order[ordinary + 1]
    if np.any(corners[first] == corners[second]):
        raise ValueError("Adjacent faces have inconsistent winding")
    extra = []
    for start, count in zip(cuts[:-1][counts != 2], counts[counts != 2]):
        group = order[start:start+count]
        edge = (low[group[0]], high[group[0]])
        if count < 2 or count % 2:
            raise ValueError(f"Open boundary along edge {edge}")
        extra.extend(_pair_junction(edge, group, corners, coords))
    if extra:
        a, b = np.array(extra, dtype=np.int64).T
        first = np.concatenate((first, a))
        second = np.concatenate((second, b))

    # The link of a geometric vertex may consist of several disjoint cycles.
    # Triangle corners on one cycle become one indexed vertex.
    row = np.concatenate((first, second))
    col = np.concatenate((next_corner[second], next_corner[first]))
    graph = coo_matrix((np.ones(len(row), dtype=np.int8), (row, col)),
                       shape=(len(corners), len(corners))).tocsr()
    n, ids = connected_components(graph, directed=False)
    old = np.empty(n, dtype=np.int64)
    old[ids] = corners
    xyz = np.array([to_world(coords[i]) for i in old], dtype=np.float64)
    faces = ids.reshape(-1, 3)
    if not extra:
        return xyz, faces

    # Two distinct edge copies can have identical indexed endpoints when both
    # sheets reconnect at the ends of their contact. Split each such edge copy.
    e0 = np.minimum(ids[first], ids[next_corner[first]])
    e1 = np.maximum(ids[first], ids[next_corner[first]])
    ordered = np.lexsort((e1, e0))
    repeated = (np.diff(e0[ordered]) == 0) & (np.diff(e1[ordered]) == 0)
    duplications = np.unique(np.r_[ordered[:-1][repeated], ordered[1:][repeated]])
    if not len(duplications):
        return xyz, faces
    middle = np.full(len(corners), -1, dtype=np.int64)
    more = []
    for e in duplications:
        a, b = first[e], second[e]
        middle[a] = middle[b] = n + len(more)
        more.append((xyz[ids[a]] + xyz[ids[next_corner[a]]]) * 0.5)
    affected = np.any(middle.reshape(-1, 3) != -1, axis=1)
    new_faces = []
    for i in np.flatnonzero(affected):
        center = n + len(more)
        more.append(xyz[faces[i]].mean(axis=0))
        ring = []
        for j in range(3*i, 3*i+3):
            ring.append(ids[j])
            if middle[j] != -1:
                ring.append(middle[j])
        for a, b in zip(ring, ring[1:]+ring[:1]):
            new_faces.append((a, b, center))
    return np.vstack((xyz, more)), np.vstack((faces[~affected], new_faces))


def stitch(xy, pieces, to_world=lambda point: point):
    """Build boundary of a union of ranked air spans over conforming CCW triangles.

    Heights at each shared projected vertex must use the same exact coordinate
    system and footprint vertex IDs. Distinct stock sheets retain their exact
    heights until the complete boundary has been indexed.
    """
    heights = defaultdict(dict)
    coords = []
    faces = []
    edges = defaultdict(list)

    def vertex(i, z):
        column = heights[i]
        if z not in column:
            column[z] = len(coords)
            coords.append((xy[i][0], xy[i][1], z))
        return column[z]

    for i, j, k, lo, hi in pieces:
        if lo == hi:
            continue
        a, b, c = (vertex(node, z) for node, z in zip((i, j, k), lo))
        faces.append((a, c, b))
        a, b, c = (vertex(node, z) for node, z in zip((i, j, k), hi))
        faces.append((a, b, c))
        for p, q, l0, l1, h0, h1 in ((i, j, lo[0], lo[1], hi[0], hi[1]),
                                       (j, k, lo[1], lo[2], hi[1], hi[2]),
                                       (k, i, lo[2], lo[0], hi[2], hi[0])):
            if l0 == h0 and l1 == h1:
                continue
            if p < q:
                edges[p, q].append(((l0, l1), (h0, h1), 1))
            else:
                edges[q, p].append(((l1, l0), (h1, h0), -1))
    if not faces:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int64)

    columns = {i: sorted(h) for i, h in heights.items()}
    for (i, j), spans in edges.items():
        if len(spans) == 2 and spans[0][:2] == spans[1][:2] and spans[0][2] == -spans[1][2]:
            continue
        events = defaultdict(int)
        for lo, hi, side in spans:
            events[lo] += side
            events[hi] -= side
        levels = sorted(events, key=lambda z: z[0] + z[1])
        inside = 0
        for lo, hi in zip(levels, levels[1:]):
            inside += events[lo]
            if not inside:
                continue
            if abs(inside) != 1 or hi[0] < lo[0] or hi[1] < lo[1]:
                raise ValueError("Air spans cross or overlap along a projected edge")
            # Each wall face must use all corner-column levels so even a
            # vertical T-junction has identical subdivision on its two sides.
            left = [heights[i][z] for z in columns[i][bisect_left(columns[i], lo[0]):bisect_right(columns[i], hi[0])]]
            right = [heights[j][z] for z in columns[j][bisect_left(columns[j], lo[1]):bisect_right(columns[j], hi[1])]]
            if inside < 0:
                left, right = right, left
            faces.extend((left[0], a, b) for a, b in zip(right, right[1:]))
            faces.extend((a, right[-1], b) for a, b in zip(left, left[1:]))
    return _normalize_boundary(coords, faces, to_world)


def write_ply(path, xyz, faces):
    header = ("ply\nformat binary_little_endian 1.0\n"
              f"element vertex {len(xyz)}\n"
              "property double x\nproperty double y\nproperty double z\n"
              f"element face {len(faces)}\n"
              "property list uchar int vertex_indices\nend_header\n")
    records = np.empty(len(faces), dtype=[("count", "u1"), ("indices", "<i4", (3,))])
    records["count"] = 3
    records["indices"] = faces
    with path.open("wb") as stream:
        stream.write(header.encode("ascii"))
        np.ascontiguousarray(xyz, dtype="<f8").tofile(stream)
        records.tofile(stream)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=True, parents=True)
    vertices, faces, scale = stock(args.part)
    setups = (("axial", (0, 0, 1)), ("side", (1, 0, 0)), ("inclined", (2, -1, 4)))
    with TemporaryDirectory(dir=args.out) as temp:
        for name, direction in setups:
            projected, to_world = pull_frame(vertices, direction, scale)
            xyz, facets = air_boundary(projected, faces, scale, to_world)
            write_ply(Path(temp) / (name + ".ply"), xyz, facets)
        for name, _ in setups:
            filename = name + ".ply"
            (Path(temp) / filename).replace(args.out / filename)


if __name__ == "__main__":
    main()
PULL_PYTHON
chmod 755 /app/out/pull.py
/opt/conda/bin/python /app/out/pull.py --part "/app/inputs/hexapod_spherical_robot/printed_parts/upper_spherical_shell.STL" --out /app/out
