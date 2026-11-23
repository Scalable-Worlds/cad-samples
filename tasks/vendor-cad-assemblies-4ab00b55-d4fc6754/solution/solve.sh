#!/bin/bash
set -euo pipefail
mkdir -p /app/out
for deliverable in /app/out/cutblank.py /app/out/blank.dxf; do
    if [[ -f "$deliverable" || -L "$deliverable" ]]; then
        unlink -- "$deliverable"
    fi
done
cat > /app/out/cutblank.py <<'CUTBLANK_PY'
#!/opt/conda/bin/python
"""Develop the picked skin of a 1.5 mm formed sheet from a binary STL into a 1:1 DXF."""
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def unit_vector(v):
    return v / np.linalg.norm(v)


def mesh_edges(faces):
    pairs = faces[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2)
    return np.unique(np.sort(pairs, axis=1), axis=0,
                     return_inverse=True, return_counts=True)


def read_sheet_stl(path):
    """Build a consistently stitched triangle mesh from the STL's facet soup."""
    record = np.dtype([('normal', '<f4', 3), ('points', '<f4', (3, 3)),
                       ('attribute', '<u2')])
    with open(path, 'rb') as stream:
        stream.seek(80)
        count = int.from_bytes(stream.read(4), 'little')
        data = np.fromfile(stream, dtype=record, count=count)
    if not count or len(data) != count:
        raise ValueError('Incomplete binary STL')
    points, indices = np.unique(data['points'].reshape(-1, 3), axis=0,
                                return_inverse=True)
    vertices = points.astype(float)
    faces = indices.reshape(-1, 3)
    _, inverse, counts = mesh_edges(faces)
    inverse = inverse.ravel()

    if np.any(counts == 1):
        # Preserve already closed trims. Match free half-edges in reverse
        # order using both endpoints, including seams with only one split end.
        open_edges = np.flatnonzero(counts[inverse] == 1)
        starts = faces.ravel()[open_edges]
        ends = np.roll(faces, -1, axis=1).ravel()[open_edges]
        midpoint = (vertices[starts] + vertices[ends]) / 2
        # Keep the stitching tolerance independent of the model translation.
        weld = 1.000001e-6
        candidates = cKDTree(midpoint).query_pairs(weld, output_type='ndarray')
        first, second = candidates.T
        mismatch = np.maximum(
            np.linalg.norm(vertices[starts[first]] - vertices[ends[second]], axis=1),
            np.linalg.norm(vertices[ends[first]] - vertices[starts[second]], axis=1))
        possible = (mismatch <= weld) & (open_edges[first] // 3 != open_edges[second] // 3)
        candidates, mismatch = candidates[possible], mismatch[possible]
        roots = np.arange(len(vertices))

        def root(index):
            while index != roots[index]:
                roots[index] = roots[roots[index]]
                index = roots[index]
            return index

        paired = np.zeros(len(open_edges), dtype=bool)
        for left, right in candidates[np.argsort(mismatch)]:
            if paired[left] or paired[right]:
                continue
            roots[root(ends[right])] = root(starts[left])
            roots[root(starts[right])] = root(ends[left])
            paired[[left, right]] = True
        if not paired.all():
            raise ValueError('The open STL seams could not be paired')
        _, group = np.unique([root(i) for i in range(len(vertices))], return_inverse=True)
        population = np.bincount(group)
        vertices = np.column_stack([np.bincount(group, weights=coord) / population
                                    for coord in vertices.T])
        faces = group[faces]
        _, inverse, counts = mesh_edges(faces)
        inverse = inverse.ravel()

    if np.any(counts != 2):
        raise ValueError('Sheet boundaries must be shared by two triangles')
    occurrences = np.argsort(inverse).reshape(-1, 2)
    neighbor = np.empty(faces.size, dtype=int)
    neighbor[occurrences] = occurrences[:, ::-1] // 3
    return vertices, faces, neighbor.reshape(faces.shape), data["points"]


def flood_facets(seed, admissible, neighbor):
    visited = {int(seed)}
    stack = [int(seed)]
    while stack:
        for face in neighbor[stack.pop()]:
            if admissible[face] and face not in visited:
                visited.add(face)
                stack.append(face)
    return np.fromiter(visited, dtype=int)


def exposed_arc_points(x, y, tolerance):
    """Find circle support knots, omitting vertices added inside STL chords."""
    pts = np.unique(np.column_stack((x, y)), axis=0)
    ordering = pts.sum(axis=1)
    first, last = int(ordering.argmin()), int(ordering.argmax())
    chord = pts[last] - pts[first]
    span = np.linalg.norm(chord)
    if span == 0:
        return pts[[first]]
    cross = chord[0]*(pts[:, 1]-pts[first, 1]) - chord[1]*(pts[:, 0]-pts[first, 0])
    if max(abs(cross)) <= tolerance * span:
        return pts[[first, last]]
    # x+y increases strictly along this quarter circle. In this order,
    # the exposed arc is the lower, counterclockwise chain of the 2-D hull.
    hull = ConvexHull(pts).vertices
    order = pts[hull].sum(axis=1)
    hull = np.roll(hull, -int(order.argmin()))
    chain = pts[hull[:int(np.argmax(pts[hull].sum(axis=1)))+1]]
    support = []
    for point in chain:
        while len(support) >= 2:
            base = support[-2]
            a, b = support[-1]-base, point-base
            if abs(a[0]*b[1]-a[1]*b[0]) > tolerance*np.linalg.norm(b):
                break
            support.pop()
        support.append(point)
    return np.array(support)


def panel_plane(stored_facets, face_ids):
    """Orthogonal least squares with respect to area on the original STL triangles."""
    triangles = stored_facets[face_ids].astype(float)
    oriented = np.cross(triangles[:, 1] - triangles[:, 0],
                        triangles[:, 2] - triangles[:, 0])
    weights = np.linalg.norm(oriented, axis=1)
    weights /= weights.sum()
    centers = triangles.mean(axis=1)
    origin = centers[0] + weights @ (centers - centers[0])
    # Equal-weight quadrature at barycentric (2/3, 1/6, 1/6) and its two
    # permutations integrates every quadratic on a triangle exactly. In
    # Cartesian space the samples are halfway from centroid to vertex.
    offsets = (centers - origin)[:, None, :] + 0.5 * (triangles - centers[:, None, :])
    offsets *= np.sqrt(weights)[:, None, None]
    # The singular vector minimizes the integrated squared distance. QR
    # keeps the final SVD small without squaring a slender panel's condition number.
    _, _, directions = np.linalg.svd(np.linalg.qr(offsets.reshape(-1, 3), mode='r'))
    normal = directions[-1]
    if normal @ oriented.sum(axis=0) < 0:
        normal = -normal
    return origin, normal


def segment(part, datum):
    v, f, nbr, stored_facets = read_sheet_stl(part)
    datum = np.asarray(datum, dtype=float)
    p = v[f]
    c = p.mean(axis=1)
    a, b = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    area_vectors = np.cross(a, b)
    area = np.linalg.norm(area_vectors, axis=1)
    if np.any(area == 0):
        raise ValueError('STL contains a degenerate triangle')
    normals = area_vectors / area[:, None]
    altitude = area / np.maximum.reduce([np.linalg.norm(a, axis=1),
                                         np.linalg.norm(b, axis=1),
                                         np.linalg.norm(a-b, axis=1)])
    ulps = np.spacing(np.abs(v).astype('f4')).astype(float)
    # The datum may project exactly onto an internal triangulation edge.
    d = datum - p[:, 0]
    s = np.einsum('ij,ij->i', np.cross(d, b), normals) / area
    t = np.einsum('ij,ij->i', np.cross(a, d), normals) / area
    dist = np.abs(np.einsum('ij,ij->i', d, normals))
    dist[(s < -1.e-5) | (t < -1.e-5) | (s+t > 1+1.e-5)] = np.inf
    root_face = int(dist.argmin())
    if not np.isfinite(dist[root_face]):
        raise ValueError('Could not locate the datum on the sheet')
    # Absolute positional allowance for vertex rounding in a float STL.
    max_ulp = ulps.max(axis=0)
    face_eps = 3.e-6 + 2 * (abs(normals) @ max_ulp)
    # Recover flat facets subdivided by STL editing or refinement. A trim
    # wall patch reaches the other stock skin even when an individual small
    # triangle on that wall happens to lie inside the bend's radius band.
    adj = nbr.ravel()
    fid = np.repeat(np.arange(len(f)), 3)
    lo, hi = fid[fid < adj], adj[fid < adj]
    delta = c[hi]-c[lo]
    flat = ((np.einsum('ij,ij->i', normals[hi], normals[lo]) > .99985)
            & (abs(np.einsum('ij,ij->i', normals[lo], delta)) < face_eps[lo])
            & (abs(np.einsum('ij,ij->i', normals[hi], delta)) < face_eps[hi]))
    graph = coo_matrix((np.ones(flat.sum(), dtype=bool), (lo[flat], hi[flat])),
                       shape=(len(f), len(f)))
    npatch, patch = connected_components(graph, directed=False)
    patch_vector = np.zeros((npatch, 3))
    np.add.at(patch_vector, patch, area_vectors)
    patch_area = np.linalg.norm(patch_vector, axis=1)
    patch_normal = patch_vector / np.maximum(patch_area[:, None], 1.e-30)
    charts, parents = [], []
    owner = np.full(len(f), -1, dtype=int)

    def append(chart, parent):
        idx = len(charts)
        owner[chart['faces']] = idx
        charts.append(chart)
        parents.append(parent)
        return idx

    def add_plane(seed, parent, hint):
        # A fan of long, thin triangles may separate two well-conditioned
        # regions of one plane. Test normal disagreement as a displacement
        # over each facet's altitude, allowing for its float32 uncertainty.
        projection_error = 2.e-6 + max_ulp @ abs(hint)
        cosine = np.clip(normals @ hint, -1., 1.)
        disagreement = np.sqrt(2 * (1 - cosine)) * altitude
        aligned = ((owner == -1) & (cosine > 0.)
                   & (disagreement <= .01 * altitude + 3 * projection_error))
        near = flood_facets(seed, aligned, nbr)
        # Facets with the largest altitudes give the most reliable normals.
        # Near-tangent chords can also lie in the coarse component. Filter
        # those out before estimating the plane used to grow the whole panel.
        best = near[np.argmax(altitude[near])]
        n = normals[best]
        good = near[(normals[near] @ n > 1 - 2.e-7)
                    & (altitude[near] >= altitude[best] * .04)]
        point, n = panel_plane(stored_facets, good)
        # On a short land, rounding error divided by the transverse span
        # can exceed the angular accuracy of the arriving quarter. Use the
        # tangent to locate such a panel and retain the measured fit on
        # well-conditioned panels.
        use_guide = False
        if parent >= 0:
            across = np.cross(charts[parent]['axis'], hint)
            uncertainty = 3.e-6 + 2 * (max_ulp @ abs(hint))
            use_guide = np.ptp((p[good] - point) @ across) * 2.e-5 <= uncertainty
            if use_guide:
                point, n = c[seed], hint
        # Use each rounded vertex's uncertainty in the plane normal direction.
        # In particular an axis-aligned panel can have an almost exact normal
        # coordinate even when the other coordinates are hundreds of mm.
        band = 2.e-6 + 1.2 * (ulps @ abs(n))
        if use_guide:
            # Include the uncertainty of the tangent-plane anchor.
            band += 2.e-6 + 1.2 * (max_ulp @ abs(n))
        residual = np.abs((v-point) @ n)
        # The orientation of a sliver has angular error proportional to its
        # vertex position uncertainty divided by the minimum altitude.
        projected_noise = 2.e-6 + max_ulp @ abs(n)
        facing = (normals @ n > .95) | ((normals @ n > 0.) & (altitude < 8*projected_noise))
        aligned = ((owner == -1) & facing
                   & (residual[f] <= band[f] + 3.e-7).all(axis=1))
        ids = flood_facets(seed, aligned, nbr)
        # The drawing normal is the area-integral least-squares normal over
        # the complete connected panel, including slivers and vertex scatter.
        point, n = panel_plane(stored_facets, ids)
        if parent >= 0 and np.ptp((p[ids] - point) @ across) * 2.e-5 <= uncertainty:
            n = hint
        return append(dict(kind='plane', faces=ids, normal=n, point=point), parent)

    def boundary(chart, include_mapped=False):
        faces = chart['faces']
        other = nbr[faces]
        adjacent = owner[other]
        candidates = (adjacent != owner[faces, None]) if include_mapped else (adjacent == -1)
        i, j = np.where(candidates)
        i = faces[i]
        va, vb = f[i, j], f[i, (j+1) % 3]
        length = np.linalg.norm(v[va]-v[vb], axis=1)
        order = np.argsort(-length)
        return zip(nbr[i,j][order], va[order], vb[order])

    def add_cylinder(seed, va, vb, parent, seams, partials, reverse=False):
        panel = charts[parent]
        n = panel['normal']
        # At a tangent seam the mesh edge is an axial generator. The two
        # signed transverse coordinates of all bend vertices must fit a
        # quarter circle tangent to the known panel, of an admissible radius.
        if normals[seed] @ n < 0.55:
            return
        origin = (v[va] + v[vb]) / 2
        u = v[vb] - v[va]
        u = unit_vector(u - n*(u @ n))
        eps = max(5.e-6, 2.5 * max_ulp @ abs(np.cross(n, u)))
        # Extend the seed through subdivisions and gaps in the same straight
        # tangency. A short seed chord is less reliable at a distant station.
        # Use the two extreme seam vertices, rather than weighting the fit
        # by how densely particular pieces of the seam were triangulated.
        # Consistent sheet winding points along the tangent with material
        # on the same side. The other edge of a short land runs in reverse.
        forward = (v[seams[:,1]] - v[seams[:,0]]) @ u > 0
        seams = seams[forward]
        seam_offsets = v[seams]
        span = np.linalg.norm(v[vb] - v[va])
        previous = None
        for _ in range(8):
            lateral = unit_vector(np.cross(n, u))
            offsets = seam_offsets - origin
            stations = offsets @ u
            allowance = eps * (1.0 + 2.0 * abs(stations) / span)
            selected = (abs(offsets @ lateral) <= allowance).all(axis=1)
            edge_ids = np.unique(seams[selected])
            if len(edge_ids) < 2:
                break
            station = (v[edge_ids] - origin) @ u
            lo, hi = edge_ids[station.argmin()], edge_ids[station.argmax()]
            if previous == (lo, hi):
                break
            previous = (lo, hi)
            chord = v[hi] - v[lo]
            origin = (v[hi] + v[lo]) / 2
            chord -= n * (chord @ n)
            span = np.linalg.norm(chord)
            u = chord / span
        origin = origin - n * ((origin-panel['point']) @ n)
        xdir = np.cross(n, u)
        step = c[seed] - v[va]
        if step @ xdir < 0:
            xdir = -xdir
        sign = 1 if step @ n < 0 else -1
        if reverse:
            sign = -sign
        # Physical radius range of this skin (inside + gauge on a convex pick).
        gauge_offset = 1.5 if sign > 0 else 0.
        rmin, rmax = 1. + gauge_offset, 2. + gauge_offset
        chord_error = .09
        x = (v-origin) @ xdir
        y = -sign * ((v-origin) @ n)
        normal_eps = 2.e-6 + 1.5 * (max_ulp @ abs(n))
        margin = max(.0001, 1.5 * eps, 2 * normal_eps)
        # A near-tangent chord can be indistinguishable from its panel.
        # If the panel contains such chords, the seed generator is later
        # around the quarter arc by up to sqrt(2*r*normal uncertainty).
        tangent_shift = np.sqrt(2 * rmax * (2*normal_eps))
        # Start with a permissive radius envelope, fit the actual tangency
        # position as well as the radius, then collect the fitted quarter.
        rr = x*x+y*y
        rough = chord_error + margin + tangent_shift
        wedge = ((x >= -rough) & (y >= -margin) & (x >= y-rough)
                 & (x <= rmax+rough) & (rr <= 2*rmax*y+8*rough))
        # A cylinder chord may itself be subdivided into many coplanar
        # facets. Their new vertices lie below the true circular radius.
        # Across a chord the geometric normal still points into the same
        # radial quadrant. End and cut walls instead face along the axis
        # or transverse to the radial field.
        nx, ny = sign*(normals @ xdir), normals @ n
        xx, yy = x[f].mean(axis=1), y[f].mean(axis=1)
        rm = (rmin+rmax)/2
        facing = ((nx >= -.04) & (ny >= -.04)
                  & (nx*xx + ny*(rm-yy) > .49*np.hypot(xx,rm-yy)))
        # Retain even tilted cylinder triangles when all vertex samples lie
        # in the original admissible analytic annulus, as may occur along
        # an oblique trim of the cylinder.
        ring = (rr >= 2*rmin*y-8*rough)
        fit_seed = facing | (ring[f].all(axis=1) & (abs(normals @ u) < .98)
                             & ((nx + ny) > .12))
        ok = (owner == -1) & wedge[f].all(axis=1) & fit_seed
        if not ok[seed]:
            return
        ids = flood_facets(seed, ok, nbr)
        verts = np.unique(f[ids])
        # Find the actual tangent generator on the far panel. Its distance
        # to the entry plane is the skin radius. This construction uses
        # both tangent planes, so vertices inserted in STL chords cannot
        # bias the radius on either skin. Walk from facets spanning the
        # bend towards a coplanar exit patch, across a generator.
        a = f[ids]
        b = np.roll(a, -1, axis=1)
        next_face = nbr[ids]
        # x,y are transverse; axial stations may differ by the full width.
        edge_x = (x[a]+x[b])/2
        edge_y = (y[a]+y[b])/2
        seam_mask = ((abs(x[a]-x[b]) <= 2*eps)
                     & (abs(y[a]-y[b]) <= 2*eps)
                     & (edge_y >= rmin-rough) & (edge_y <= rmax+rough)
                     & (abs(edge_x-edge_y) <= 2*rough)
                     & ((nx[next_face] > .999)
                        | ((nx[next_face] > 0.) & (altitude[next_face]
                           * np.sqrt(np.maximum(0., 2*(1-nx[next_face]))) <= 3*face_eps[next_face])))
                     & (np.ptp(x[f[next_face]], axis=-1) <= 3*eps)
                     & (yy[next_face] > edge_y)
                     & (yy[ids, None] < edge_y))
        rows, cols = np.where(seam_mask)
        if len(rows):
            # An almost tangent last chord can project like the exit plane
            # after STL rounding. The exit patch has the orthogonal normal.
            # Estimate it from the full area-weighted patch, not one sliver.
            target_patch = patch[next_face[rows,cols]]
            # A fine mesh can put the first few chords inside the plane's
            # rounding band. Recover the exit station from the circular
            # support between the two tangent planes. If z is the distance
            # to the exit plane and y the distance to the entry plane, a
            # point on the quarter satisfies r = y + z + sqrt(2*y*z).
            exit_x, exit_y = edge_x[rows,cols], edge_y[rows,cols]
            mismatch = abs(exit_x - exit_y)
            if np.ptp(exit_y) > .65*(eps + normal_eps):
                samples = exposed_arc_points(x[verts], y[verts], max(2.e-6, .7*eps))
                z = np.median(exit_x) - samples[:,0]
                rise = samples[:,1]
                middle = (z > .03*rmin) & (rise > .03*rmin)
                if middle.any():
                    z, rise = z[middle], rise[middle]
                    r_support = np.median(z + rise + np.sqrt(2*z*rise))
                    mismatch = abs(exit_y - r_support)
            feasible = mismatch <= mismatch.min() + .65*(eps+normal_eps)
            direction_error = np.linalg.norm(np.cross(patch_normal[target_patch], xdir), axis=1)
            feasible &= direction_error <= np.min(direction_error[feasible]) + 2.e-5
            weight = np.where(feasible, patch_area[target_patch], 0.)
            primary = target_patch[np.argmax(weight)]
            lengths = np.linalg.norm(v[a[rows,cols]]-v[b[rows,cols]], axis=1)
            k = int(np.argmax(np.where(feasible & (target_patch == primary), lengths, -1.)))
            same_generator = ((abs(exit_x-exit_x[k]) <= 2*eps)
                              & (abs(exit_y-exit_y[k]) <= 2*eps))
            near_rows, near_cols = rows[same_generator], cols[same_generator]
            ends = np.unique(np.r_[a[near_rows,near_cols], b[near_rows,near_cols]])
            stations = (v[ends]-origin) @ u
            ends = ends[[stations.argmin(), stations.argmax()]]
            r = y[ends].mean()
            h = x[ends].mean() - r
        else:
            if not partials:
                return
            # A relief can isolate a tongue on a bend without leaving it
            # an exit tangent. Use the complete sibling's geometry when
            # its entry line and circular support agree with this tongue.
            samples = exposed_arc_points(x[verts], y[verts], max(2.e-6, eps*.7))
            xx_sample, yy_sample = samples.T
            fits = []
            for sibling in charts:
                if sibling['kind'] != 'cylinder' or sibling['sign'] != sign:
                    continue
                if abs(sibling['axis'] @ u) < 1 - 1.e-8:
                    continue
                r0 = sibling['radius']
                d = sibling['center'] - origin
                # Compare the two axes at this tongue's station.
                d = d - (d @ sibling['axis']) * sibling['axis']
                h0 = d @ xdir
                if abs(d @ n + sign*r0) > 2*margin or abs(h0) > tangent_shift + 2*margin:
                    continue
                error = max(abs(np.hypot(xx_sample-h0, yy_sample-r0) - r0))
                if error <= 2*margin:
                    fits.append((error, h0, r0))
            if fits:
                _, h, r = min(fits)
            else:
                outer = yy_sample > max(eps, .12 * yy_sample.max())
                if not outer.any():
                    return
                # With entry tangent fixed, x*x+y*y = 2*r*y.
                estimates = (xx_sample[outer]**2+yy_sample[outer]**2) / (2*yy_sample[outer])
                h, r = 0., float(estimates.max())
        if not rmin - 2*margin <= r <= rmax + 2*margin:
            return
        radial = np.hypot(x-h, y-r)
        quarter = (x-h >= -margin) & (y <= r+margin)
        on_circle = (abs(radial-r) <= margin)[f].all(axis=1)
        chord_band = (radial <= r+margin) & (radial >= r-chord_error-margin)
        close_normal = (nx[:,None]*(x[f]-h) + ny[:,None]*(r-y[f]) > .61*radial[f])
        in_quarter = quarter[f].all(axis=1)
        in_band = chord_band[f].all(axis=1) & in_quarter
        escape = np.zeros(npatch, dtype=bool)
        escape[patch[~in_band]] = True
        on_chord = in_band & close_normal.all(axis=1) & ~escape[patch]
        ok = (owner == -1) & in_quarter & (on_circle | on_chord)
        if not ok[seed]:
            return
        if len(rows):
            # The far tangent is an observed interface. Stop the cylinder
            # there even if a very short land fits its positional error band.
            ok[next_face[near_rows,near_cols]] = False
        ids = flood_facets(seed, ok, nbr)
        axis = unit_vector(np.cross(sign*n, xdir))
        center = origin + h*xdir - sign*r*n
        append(dict(kind='cylinder', faces=ids, axis=axis, center=center,
                    radius=r, sign=sign, start_normal=n), parent)

    # A narrow triangle containing the datum can have a poorly resolved
    # normal. Find its best-conditioned neighbour on the planar face.
    # The datum clearance and trim chord bound leave this 0.35 mm ball
    # inside the face. Traverse only edges which actually meet the ball.
    local = {root_face}
    todo = [root_face]
    while todo:
        face = todo.pop()
        start = p[face]
        edge = np.roll(start, -1, axis=0)-start
        along = np.einsum('ij,ij->i', datum-start, edge) / np.einsum('ij,ij->i', edge, edge)
        nearest = start + np.clip(along, 0., 1.)[:,None]*edge
        for adjacent in nbr[face, np.linalg.norm(nearest-datum, axis=1) <= .35]:
            if adjacent not in local:
                local.add(adjacent)
                todo.append(adjacent)
    root_face = max(local, key=altitude.__getitem__)
    add_plane(root_face, -1, normals[root_face])
    # Discover connected panels through full quarters first. This gives
    # even a tongue reached on an early panel access to a sibling cylinder
    # recognized later in the traversal.
    for partials in (False, True):
        for idx, chart in enumerate(charts):
            if chart['kind'] == 'plane':
                candidates = list(boundary(chart, include_mapped=True))
                seams = np.array([(a,b) for i,a,b in candidates if normals[i] @ chart['normal'] > .55], dtype=int)
                for seed, va, vb in candidates:
                    if owner[seed] == -1:
                        add_cylinder(seed, va, vb, idx, seams, partials)
                        # On an almost tangent first chord float rounding can
                        # reverse the transverse displacement of its centroid.
                        if (owner[seed] == -1 and normals[seed] @ chart['normal'] > .55
                                and abs((c[seed]-v[va]) @ chart['normal']) <= face_eps[seed]):
                            add_cylinder(seed, va, vb, idx, seams, partials, reverse=True)
            else:
                n0, axis, ctr = chart['start_normal'], chart['axis'], chart['center']
                n1 = np.cross(axis, n0)
                # A fork or slot can leave several disjoint panel patches at
                # either tangent. Retain their actual edge adjacency.
                for seed, va, vb in boundary(chart):
                    if owner[seed] != -1:
                        continue
                    for n in (n1, n0):
                        dot = normals[seed] @ n
                        angular = np.linalg.norm(np.cross(normals[seed], n))
                        if dot <= 0 or (dot < .9998 and angular*altitude[seed] > 3*face_eps[seed]):
                            continue
                        tangent_pt = ctr + chart['sign'] * chart['radius'] * n
                        lateral = np.cross(axis, n)
                        offsets = v[[va,vb]]-tangent_pt
                        if np.max(abs(offsets @ lateral)) < .01 and np.max(abs(offsets @ n)) < .01:
                            add_plane(seed, idx, n)
                            break
    return v, f, charts, parents, 0


# The maps use row vectors for points. A plane basis is a pair of rows
# taking a spatial displacement to the two drawing coordinates.
class _PlaneMap:
    def __init__(self, base, origin, basis):
        self.base = base
        self.origin = origin
        self.basis = basis

    def __call__(self, points):
        return self.origin + (points - self.base) @ self.basis.T


class _BendMap:
    def __init__(self, chart, plane, normal, vertices):
        self.center = chart["center"]
        self.axis = unit_vector(chart["axis"])
        self.sign = chart["sign"]
        self.radius = chart["radius"]
        self.neutral_radius = self.radius - self.sign * 0.75
        self.radial = self.radial_at_plane(normal)
        self.tangent = np.cross(self.axis, self.radial)

        # The plane touches the chosen skin here. Shifting both the plane
        # and this tangent generator by -t/2 along n leaves its xy unchanged.
        contact = self.center + self.radius * self.radial
        self.origin = plane(contact)
        self.axial_xy = unit_vector(plane.basis @ self.axis)
        turn = plane.basis @ self.tangent
        self.arc_xy = unit_vector(turn - np.dot(turn, self.axial_xy) * self.axial_xy)

        # With a slotted bend the angular domain need not be populated at
        # every axial station. All of it lies on one side of the tangent.
        angles, _ = self.parameters(vertices)
        self.sweep = np.copysign(np.pi / 2, angles[np.argmax(abs(angles))])

    def radial_at_plane(self, normal):
        return self.sign * unit_vector(normal - np.dot(normal, self.axis) * self.axis)

    def parameters(self, points):
        d = points - self.center
        return np.arctan2(d @ self.tangent, d @ self.radial), d @ self.axis

    def xy(self, theta, station):
        return (self.origin + np.asarray(station)[..., None] * self.axial_xy
                + (self.neutral_radius * np.asarray(theta))[..., None] * self.arc_xy)

    def __call__(self, points):
        return self.xy(*self.parameters(points))

    def attached_plane(self, normal):
        # A fork may reconnect to the entry panel on another interval.
        # The fitted face normal identifies the tangent; its measurement
        # noise must not change the specified quarter-turn angle.
        entry = self.radial
        exit_radial = np.sign(self.sweep) * self.tangent
        if self.sign * (normal @ (exit_radial - entry)) > 0:
            angle, radial = self.sweep, exit_radial
        else:
            angle, radial = 0., entry
        tangent = np.cross(self.axis, radial)
        base = self.center + self.radius * radial
        basis = np.outer(self.axial_xy, self.axis) + np.outer(self.arc_xy, tangent)
        return _PlaneMap(base, self.xy(angle, 0), basis)


def _develop_maps(vertices, triangles, charts, parents, root, datum):
    n = unit_vector(charts[root]["normal"])
    reference = np.eye(3)[np.argmin(np.round(abs(n), 6))]
    x = unit_vector(reference - np.dot(reference, n) * n)
    maps = [None] * len(charts)
    maps[root] = _PlaneMap(datum, np.zeros(2), np.array([x, np.cross(n, x)]))

    # Charts are numbered in discovery order. Each parent already has a
    # drawing map when its child's tangent generator is developed.
    for child, parent in enumerate(parents):
        if child == root:
            continue
        chart = charts[child]
        if chart["kind"] == "plane":
            maps[child] = maps[parent].attached_plane(chart["normal"])
        else:
            ids = np.unique(triangles[chart["faces"]])
            maps[child] = _BendMap(chart, maps[parent], charts[parent]["normal"],
                                   vertices[ids])
    return maps


def _cut_loops(skin, xy):
    # Retain only free half-edges of the chosen continuous skin. This
    # cancels every triangulation edge and every internal tangent seam.
    edges = np.vstack((skin[:, [0, 1]], skin[:, [1, 2]], skin[:, [2, 0]]))
    _, first, count = np.unique(np.sort(edges, axis=1), axis=0,
                                return_index=True, return_counts=True)
    boundary = edges[first[count == 1]]
    following = dict(boundary.tolist())
    if len(following) != len(boundary) or set(following) != set(following.values()):
        raise ValueError("Skin trim edges do not form disjoint closed contours")

    loops = []
    while following:
        start = v = next(iter(following))
        ids = []
        while v in following:
            ids.append(v)
            v = following.pop(v)
        if v != start or len(ids) < 3:
            raise ValueError("Open or degenerate contour on the selected stock skin")
        points = xy[ids]
        keep = np.any(points != np.roll(points, 1, axis=0), axis=1)
        loops.append(points[keep])
    return loops


def _brake_lines(bend, vertices, triangles):
    """Section the mesh on the half-angle radial plane and union axial spans."""
    ids, local = np.unique(triangles.reshape(-1), return_inverse=True)
    angle = bend.sweep / 2
    transverse = np.cos(angle) * bend.tangent - np.sin(angle) * bend.radial
    offsets = vertices[ids] - bend.center
    station = offsets @ bend.axis

    # An oblique trim crosses the radial plane on a physical STL chord.
    # Linear interpolation of its cylindrical angle would bias the station.
    delta = offsets @ transverse
    edges, inverse, _ = mesh_edges(local.reshape(-1, 3))
    # In the half-angle strip trim branches cross transversely. Section
    # the original STL chords directly, without moving nearby vertices.

    # Compute a shared edge's intersection exactly once. Identical interval
    # ends then meet without a tolerance that could swallow a narrow relief.
    a, b = edges.T
    low = np.full(len(edges), np.inf)
    high = np.full(len(edges), -np.inf)
    crossings = (delta[a] * delta[b]) < 0
    u, v = a[crossings], b[crossings]
    fraction = delta[u] / (delta[u] - delta[v])
    low[crossings] = high[crossings] = station[u] + fraction * (station[v] - station[u])
    for endpoint in (a, b):
        on_section = delta[endpoint] == 0
        value = station[endpoint[on_section]]
        low[on_section] = np.minimum(low[on_section], value)
        high[on_section] = np.maximum(high[on_section], value)
    local_edges = inverse.reshape(-1, 3)
    starts = low[local_edges].min(axis=1)
    ends = high[local_edges].max(axis=1)
    merged = []
    for index in np.argsort(starts):
        if starts[index] >= ends[index]:
            continue
        if merged and starts[index] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], ends[index])
        else:
            merged.append([starts[index], ends[index]])
    layer = "OUTER" if bend.sign > 0 else "INNER"
    return [(layer, bend.xy(angle, low), bend.xy(angle, high)) for low, high in merged]


def develop(part, datum):
    datum = np.asarray(datum, dtype=float)
    vertices, triangles, charts, parents, root = segment(part, datum)
    maps = _develop_maps(vertices, triangles, charts, parents, root, datum)
    xy_sum = np.zeros((len(vertices), 2))
    visits = np.zeros(len(vertices), dtype=int)
    selected = []
    marks = []
    for chart, mapping in zip(charts, maps):
        facets = triangles[chart["faces"]]
        selected.append(facets)
        ids = np.unique(facets)
        xy_sum[ids] += mapping(vertices[ids])
        visits[ids] += 1
        if chart["kind"] == "cylinder":
            marks.extend(_brake_lines(mapping, vertices, facets))
    used = visits > 0
    xy_sum[used] /= visits[used, None]
    return _cut_loops(np.concatenate(selected), xy_sum), marks


def main():
    import argparse
    from pathlib import Path
    from tempfile import TemporaryDirectory
    import ezdxf

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", type=Path, required=True, help="binary STL in mm")
    parser.add_argument("--datum", type=float, nargs=3, required=True,
                        metavar=("X", "Y", "Z"), help="point on a planar stock face")
    parser.add_argument("--out", type=Path, required=True, help="drawing directory")
    args = parser.parse_args()
    loops, marks = develop(args.part, args.datum)
    doc = ezdxf.new("R2010", units=ezdxf.units.MM)
    doc.header["$MEASUREMENT"] = 1
    model = doc.modelspace()
    for name, color in (("CUT", 7), ("INNER", 3), ("OUTER", 1)):
        doc.layers.add(name, color=color)
    for loop in loops:
        model.add_lwpolyline(loop, format="xy", close=True, dxfattribs={"layer": "CUT"})
    for layer, a, b in marks:
        model.add_line(a, b, dxfattribs={"layer": layer})
    args.out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.out) as stage:
        path = Path(stage) / "blank.dxf"
        doc.saveas(path)
        path.replace(args.out / "blank.dxf")


if __name__ == "__main__":
    main()
CUTBLANK_PY
chmod 0755 /app/out/cutblank.py
/opt/conda/bin/python /app/out/cutblank.py --part /app/inputs/base.STL --datum 161.690406 1.5 156.623074 --out /app/out
