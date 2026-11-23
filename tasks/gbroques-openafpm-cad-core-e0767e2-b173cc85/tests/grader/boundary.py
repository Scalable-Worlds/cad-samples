"""Enclosures of full trimmed material boundaries in their installed frame."""
from itertools import product
import math

import FreeCAD as App
import Part
import numpy as np

TOL = .05


def _rigid(a, b):
    pa, pb = a.Placement, b.Placement
    a, b = a.cleaned(), b.cleaned()
    a.Placement = App.Placement()
    b.Placement = App.Placement()
    if a.exportBrepToString() != b.exportBrepToString():
        return math.inf
    box = a.BoundBox
    return max((pa.multVec(App.Vector(*p)) - pb.multVec(App.Vector(*p))).Length
               for p in product((box.XMin, box.XMax), (box.YMin, box.YMax),
                                (box.ZMin, box.ZMax)))


def _trim_arc(face, edge):
    curve, lo, hi = face.curveOnSurface(edge)
    if isinstance(curve, (Part.Geom2d.Line2d,Part.Geom2d.Line2dSegment)):
        return [(np.array([[curve.value(t).x, curve.value(t).y] for t in (lo, hi)]), np.ones(2))]
    exact = (Part.Geom2d.Circle2d, Part.Geom2d.Ellipse2d, Part.Geom2d.Hyperbola2d,
             Part.Geom2d.Parabola2d, Part.Geom2d.BSplineCurve2d, Part.Geom2d.BezierCurve2d)
    if not isinstance(curve, exact):
        return []
    arcs = curve.toBSpline(lo, hi).toBezier()
    return [(np.asarray([(p.x, p.y) for p in arc.getPoles()]), np.array(arc.getWeights()))
            for arc in arcs]


def _trim_bound(f, g, transform=None, match=False):
    """Bound both trimmed-domain differences by homotopy of corresponding UV loops.

    Each UV point whose inside/outside parity changes belongs to a strip swept
    by a pair of trimming arcs. Their common rational Bezier basis supplies an
    enclosure over the whole arc. Sign or magnitude of winding is immaterial.
    """
    ua = np.zeros(2)
    if len(f.Wires) != len(g.Wires):
        return None
    candidates = list(g.Wires)
    for wire in f.Wires:
        idx = min(range(len(candidates)), key=lambda k: (
            (candidates[k].CenterOfMass-wire.CenterOfMass).Length
            + abs(candidates[k].Length-wire.Length)/math.tau)) if match else 0
        other = candidates.pop(idx)
        left, right = wire.OrderedEdges, other.OrderedEdges
        if len(left) != len(right):
            return None
        reverse_loop = False
        if match and len(right) > 1:
            offset = min(range(len(right)), key=lambda k: (left[0].CenterOfMass-right[k].CenterOfMass).Length)
            forward = right[offset:] + right[:offset]
            backward = forward[:1] + forward[:0:-1]
            def cost(sequence):
                return sum((e.CenterOfMass-h.CenterOfMass).Length for e,h in zip(left, sequence))
            reverse_loop = cost(backward) < cost(forward)
            right = backward if reverse_loop else forward
        for e, h in zip(left, right):
            try:
                a, b = _trim_arc(f, e), _trim_arc(g, h)
            except (TypeError, ValueError, RuntimeError, Part.OCCError):
                # Some legal seams have no stored pcurve on this use of the
                # edge. The spatial enclosure below works on the trimmed face.
                return None
            if not a or len(a) != len(b):
                return None
            if ((e.Orientation == "Reversed") != (h.Orientation == "Reversed")) != reverse_loop:
                b = [(p[::-1], w[::-1]) for p, w in reversed(b)]
            for (p, w), (q, v) in zip(a, b):
                if min(w) <= 0 or min(v) <= 0:
                    return None
                while len(p) != len(q):
                    def elevate(points, weights):
                        net = np.column_stack((points*weights[:,None], weights))
                        blend = np.arange(1,len(net))[:,None] / len(net)
                        high = np.vstack((net[0], blend*net[:-1]+(1-blend)*net[1:], net[-1]))
                        return high[:,:2]/high[:,2,None], high[:,2]
                    if len(p) < len(q):
                        p,w = elevate(p,w)
                    else:
                        q,v = elevate(q,v)
                if transform is not None:
                    matrix, offset = transform
                    q = q @ matrix + offset
                # After translation about q[0], denominators stay above
                # their smallest weight. Normalization is homogeneous.
                w, v = w / w.max(), v / v.max()
                delta = np.max(abs(p-q), axis=0)
                delta += np.ptp(q, axis=0) * np.max(abs(w-v)) / min(w)
                ua = np.maximum(ua, delta)
    return ua


def _untrim(surface):
    while isinstance(surface, Part.RectangularTrimmedSurface):
        surface = surface.BasisSurface
    return surface


def _elementary_basis(face):
    """Identify planar and circular UV families by exact generating geometry."""
    actual = face.Surface
    basis = _untrim(actual)
    while isinstance(basis, Part.OffsetSurface):
        basis = _untrim(basis.BasisSurface)
    if isinstance(basis, Part.SurfaceOfExtrusion):
        curve = basis.BasisCurve
        if isinstance(curve, Part.Line):
            normal = curve.Direction.cross(basis.Direction)
            if normal.Length:
                return Part.Plane(face.valueAt(0, 0), normal)
        while isinstance(curve, Part.OffsetCurve):
            direction = curve.OffsetDirection
            curve = curve.BasisCurve
            if direction.cross(basis.Direction).Length != 0:
                return actual
        if isinstance(curve, Part.ArcOfCircle):
            curve = curve.Circle
        if isinstance(curve, Part.Circle) and curve.Axis.cross(basis.Direction).Length == 0:
            centre = (face.valueAt(0,0) + face.valueAt(math.pi,0)) * .5
            result = Part.Cylinder()
            result.Center, result.Axis = centre, basis.Direction
            result.Radius = (face.valueAt(0,0)-centre).Length
            return result
    if isinstance(basis, (Part.Plane, Part.Cylinder, Part.Cone)):
        return basis
    return actual


def _circular_region_bound(f, g):
    """Match complete planar circular trims independently of seam and wire order."""
    def rings(face):
        plane = face.Surface
        if not isinstance(plane, Part.Plane):
            return None
        result = []
        drift = 0.
        for wire in face.Wires:
            edges = wire.Edges
            if len(edges) != 1 or not isinstance(edges[0].Curve, Part.Circle):
                return None
            edge = edges[0]
            if not edge.isClosed() or not math.isclose(
                    edge.LastParameter-edge.FirstParameter, math.tau,
                    rel_tol=32*math.ulp(1.)):
                return None
            try:
                trim, lo, hi = face.curveOnSurface(edge)
            except (TypeError, ValueError, RuntimeError, Part.OCCError):
                return None
            if not isinstance(trim, Part.Geom2d.Circle2d) or not math.isclose(
                    abs(hi-lo), math.tau, rel_tol=32*math.ulp(1.)):
                return None
            uv = trim.Location
            rim = Part.Circle(face.valueAt(uv.x, uv.y), plane.Axis, trim.Radius)
            curve = edge.Curve
            # Enclose the actual 3-D edge as well as the face's own trimming
            # pcurve. They can carry different native representation tolerances.
            turn = min((rim.Axis-curve.Axis).Length, (rim.Axis+curve.Axis).Length)
            drift = max(drift, (rim.Center-curve.Center).Length + abs(rim.Radius-curve.Radius)
                        + curve.Radius*turn + rim.Radius*abs(abs(hi-lo)-math.tau))
            result.append((wire, rim))
        if not result or not face.isValid():
            return None
        outer = face.OuterWire
        return [c for _, c in sorted(result, key=lambda item: not item[0].isSame(outer))], drift

    left = rings(f)
    if left is None:
        return math.inf
    right = rings(g)
    if right is None:
        return math.inf
    left, df = left
    right, dg = right
    if len(left) != len(right):
        return math.inf
    a, b = left[0], right[0]
    normal = b.Axis if a.Axis.dot(b.Axis) >= 0 else -b.Axis
    rotate = App.Rotation(a.Axis, normal)
    # Map the first plane rigidly onto the second at the outer rim centre.
    # Every point is inside the outer disk, so this bounds its displacement.
    movement = (a.Center-b.Center).Length + a.Radius*(a.Axis-normal).Length
    if movement > TOL:
        return math.inf
    error = abs(a.Radius-b.Radius)
    holes = right[1:]
    for circle in left[1:]:
        centre = b.Center + rotate.multVec(circle.Center-a.Center)
        def rim_distance(other):
            return (centre-other.Center).Length + abs(circle.Radius-other.Radius)
        index = min(range(len(holes)), key=lambda i: rim_distance(holes[i]))
        error = max(error, rim_distance(holes.pop(index)))
    # During a linear homotopy of centres and radii, the only planar points
    # that change inclusion cross a moving rim. Each is within 'error' of
    # that rim at either end, including when the two holes have unequal radii.
    return movement + error + df + dg


def _complete_band(face):
    """Enclose the two full circular ends of one trimmed cylindrical band."""
    surface = face.Surface
    if not isinstance(surface, Part.Cylinder):
        return None
    ends = [e for e in face.Edges if not e.isSeam(face)]
    if len(ends) != 2 or not all(e.isClosed() and isinstance(e.Curve, Part.Circle) for e in ends):
        return None
    centers = []
    error = 0.
    for edge in ends:
        try:
            curve, lo, hi = face.curveOnSurface(edge)
        except (TypeError, ValueError, RuntimeError, Part.OCCError):
            return None
        if not isinstance(curve, (Part.Geom2d.Line2d, Part.Geom2d.Line2dSegment)):
            return None
        a, b = curve.value(lo), curve.value(hi)
        span = abs(b.x-a.x)
        if not math.isclose(span, math.tau, rel_tol=32*math.ulp(1.)):
            return None
        v = (a.y+b.y)/2
        centre = surface.Center + v*surface.Axis
        centers.append(centre)
        rim = edge.Curve
        tilt = min((surface.Axis-rim.Axis).Length, (surface.Axis+rim.Axis).Length)
        error = max(error, abs(b.y-a.y)/2 + surface.Radius*abs(span-math.tau),
                    (centre-rim.Center).Length + abs(surface.Radius-rim.Radius) + rim.Radius*tilt)
    if not face.isValid():
        return None
    return centers, error


def _band_bound(f, g):
    first = _complete_band(f)
    if first is None:
        return math.inf
    second = _complete_band(g)
    if second is None:
        return math.inf
    a, da = first
    b, db = second
    # Match azimuth after aligning the axis directions and interpolate between
    # corresponding cap centres. Both contributions bound the entire band.
    offset = min(max((p-q).Length for p,q in zip(a,ends)) for ends in (b,b[::-1]))
    x, y = f.Surface, g.Surface
    tilt = min((x.Axis-y.Axis).Length, (x.Axis+y.Axis).Length)
    return da + db + offset + abs(x.Radius-y.Radius) + min(x.Radius,y.Radius)*tilt


def _chart_bound(f, g):
    """Bound distance on the union of two corresponding analytic UV charts."""
    sf, sg = _elementary_basis(f), _elementary_basis(g)
    if type(sf) != type(sg) or not isinstance(sf, (Part.Plane, Part.Cylinder, Part.Cone)):
        return math.inf
    if isinstance(sf, Part.Plane):
        circular = _circular_region_bound(f, g)
        if circular <= TOL:
            return circular
    elif isinstance(sf, Part.Cylinder):
        cylinder = _band_bound(f, g)
        if cylinder <= TOL:
            return cylinder
    a, b = f.ParameterRange, g.ParameterRange
    vrange = (min(a[2], b[2]), max(a[3], b[3]))
    if isinstance(sf, Part.Plane):
        origin, normal = f.valueAt(0, 0), sf.Axis
        d = max(abs((g.valueAt(u, v)-origin).dot(normal))
                for u,v in product((b[0], b[1]), (b[2], b[3])))
        if d > TOL:
            return math.inf
        du, dv = f.valueAt(1,0)-origin, f.valueAt(0,1)-origin
        gram = np.array([[du.dot(du), du.dot(dv)], [du.dot(dv), dv.dot(dv)]])
        def project(point):
            return np.linalg.solve(gram, np.array([(point-origin).dot(du), (point-origin).dot(dv)]))
        origin_uv = project(g.valueAt(0,0))
        transform = np.array([project(g.valueAt(1,0))-origin_uv,
                              project(g.valueAt(0,1))-origin_uv]), origin_uv
        # Compare in the common orthogonal projection. The height interval
        # on g bounds separation from its planar image in either direction.
        uv = _trim_bound(f, g, transform)
        if uv is None or d + sum(uv) > TOL:
            uv = _trim_bound(f, g, transform, match=True)
        if uv is None:
            return math.inf
        in_plane = max((du*float(uv[0]) + dv*float(sign*uv[1])).Length for sign in (-1, 1))
        return math.hypot(d, in_plane)
    else:
        # C(v) + X(v) cos u + Y(v) sin u. Each coefficient is affine in v.
        d = du = 0.
        for v in vrange:
            def basis(face):
                p = face.valueAt(0, v)
                centre = (p + face.valueAt(math.pi, v)) * .5
                return centre, p-centre, face.valueAt(math.pi/2, v)-centre
            left, right = basis(f), basis(g)
            d = max(d, sum((x-y).Length for x, y in zip(left, right)))
            du = max(du, left[1].Length, right[1].Length)
    if d > TOL:
        return math.inf
    gap = _trim_bound(f, g)
    if gap is None or d+du*gap[0]+gap[1] > TOL:
        gap = _trim_bound(f, g, match=True)
    if gap is None:
        return math.inf
    # Unit-speed v generatrices; |dS/du| is bounded by the largest radius.
    return d + du * gap[0] + gap[1]


def _unit_chart(face):
    if not isinstance(face.Surface, (Part.Plane, Part.Cylinder, Part.Cone, Part.Sphere,
                                     Part.Toroid, Part.BezierSurface, Part.BSplineSurface)):
        return None
    nurbs = face.toNurbs().Faces[0]
    s = nurbs.Surface
    a, b, c, d = nurbs.ParameterRange
    s.segment(a, b, c, d)
    s.scaleKnotsToBounds(0., 1., 0., 1.)
    return nurbs, s, np.array([a,c]), np.array([b-a,d-c])


def _net(surface):
    pw = np.asarray(surface.getPolesAndWeights())
    if pw.ndim != 3 or pw.shape[2] != 4:
        raise ValueError("Degenerate control net")
    return pw[:,:,:3], pw[:,:,3]


def _derivative_limits(surface):
    p,w = _net(surface)
    centre = p[0,0]
    h = (p-centre) * w[:,:,None]
    diameter = np.max(np.linalg.norm(p-centre, axis=2))
    out = []
    for axis,degree,knots in ((0,surface.UDegree,surface.UKnotSequence),
                              (1,surface.VDegree,surface.VKnotSequence)):
        # Derivatives of rational B-splines from derivatives of homogeneous
        # control points and the lower bound on the denominator.
        ks = np.array(knots)
        n = p.shape[axis]
        delta = ks[degree+1:degree+n] - ks[1:n]
        coeff = degree / delta
        dh = np.linalg.norm(np.diff(h,axis=axis),axis=2)
        dw = abs(np.diff(w,axis=axis))
        if axis == 0:
            coeff = coeff[:,None]
        out.append(np.max((dh + diameter*dw)*coeff) / w.min())
    return np.asarray(out)


def _rational_bound(f, g):
    """Match exact rational supports using knot insertion and degree elevation."""
    charts = _unit_chart(f), _unit_chart(g)
    if None in charts:
        return math.inf
    f,s,fo,fs = charts[0]
    g,t,go,gs = charts[1]
    # A chart can exchange axes and traverse either parameter backwards.
    corners = list(product((0.,1.),repeat=2))
    points = [s.value(*p) for p in corners]
    options = []
    for perm in ((0,1),(1,0)):
        for signs in product((1,-1),repeat=2):
            mat = np.eye(2)[:,perm] * signs
            offset = (1-np.array(signs)) / 2
            mapped = [np.array(p) @ mat + offset for p in corners]
            dist = max((p-t.value(*q)).Length for p,q in zip(points,mapped))
            if dist <= TOL:
                options.append((dist,mat,offset))
    for _,mat,offset in sorted(options,key=lambda opt:opt[0]):
        # q = p @ mat + offset maps unit parameters of s into t.
        inverse = mat.T
        uv_map = (inverse / gs[:,None] * fs,
                  fo - (go / gs + offset) @ inverse * fs)
        gap = _trim_bound(f,g,uv_map,match=True)
        if gap is None:
            continue
        trim = np.dot(gap/fs, _derivative_limits(s))
        if trim > TOL:
            continue
        degs = np.array([s.UDegree,s.VDegree])
        degt = np.array([t.UDegree,t.VDegree])
        maximum = np.maximum(degs,np.abs(mat) @ degt).astype(int)
        S,T = s.copy(), t.copy()
        S.increaseDegree(*map(int,maximum))
        T.increaseDegree(*map(int,maximum @ abs(mat)))
        knots = []
        for axis,ks in enumerate((S.getUKnots(),S.getVKnots())):
            ta = int(np.argmax(abs(mat[axis])))
            sign = mat[axis,ta]
            kt = T.getUKnots() if ta == 0 else T.getVKnots()
            vals = set(ks) | {float((x-offset[ta])*sign) for x in kt}
            knots.append(sorted(vals))
        worst = 0.
        for i in range(len(knots[0])-1):
            for j in range(len(knots[1])-1):
                lo = np.array([knots[0][i],knots[1][j]])
                hi = np.array([knots[0][i+1],knots[1][j+1]])
                x,y = lo @ mat + offset, hi @ mat + offset
                ta,tb = np.minimum(x,y),np.maximum(x,y)
                sa,tg = S.copy(),T.copy()
                sa.segment(lo[0],hi[0],lo[1],hi[1])
                tg.segment(ta[0],tb[0],ta[1],tb[1])
                # Each clipped patch is a tensor product Bezier patch.
                if max(sa.NbUKnots,sa.NbVKnots,tg.NbUKnots,tg.NbVKnots)!=2:
                    return math.inf
                p,w = _net(sa)
                q,v = _net(tg)
                if mat[0,0] == 0:
                    q,v = q.swapaxes(0,1),v.T
                q = q[::int(sum(mat[0])), ::int(sum(mat[1]))]
                v = v[::int(sum(mat[0])), ::int(sum(mat[1]))]
                if p.shape != q.shape or min(w.min(),v.min())<=0:
                    return math.inf
                w,v = w/w.max(), v/v.max()
                spatial = np.max(np.linalg.norm(p-q,axis=2))
                diam = np.linalg.norm(np.ptp(q,axis=(0,1)))
                spatial += diam * np.max(abs(w-v))/w.min()
                worst = max(worst,spatial)
            if worst+trim > TOL:
                break
        if worst+trim <= TOL:
            return worst+trim
    return math.inf


def _round_chart(face):
    """Construct a complete rotationally invariant band, disk or annulus."""
    surface = face.Surface
    if isinstance(surface, Part.Cylinder):
        _, _, v, V = face.ParameterRange
        origin = surface.Center + v * surface.Axis
        ideal = Part.makeCylinder(surface.Radius, V-v, origin, surface.Axis)
        return ideal.Faces[0], origin, surface.Axis
    if isinstance(surface, Part.Plane):
        # A circular edge on its own says nothing about the outer trim.
        # Require complete centered circular loops on this plane. The later
        # trim homotopy checks their actual orientation and parameterization.
        if not face.isValid() or len(face.Wires) not in (1, 2) or any(not w.isClosed() for w in face.Wires):
            return None
        loops = []
        for wire in face.Wires:
            if not wire.Edges or not all(isinstance(e.Curve, Part.Circle) for e in wire.Edges):
                return None
            circle = wire.Edges[0].Curve
            if any(e.Curve.Center != circle.Center or e.Curve.Radius != circle.Radius
                   for e in wire.Edges):
                return None
            if not math.isclose(sum(e.Length for e in wire.Edges), math.tau*circle.Radius,
                                rel_tol=32*math.ulp(1.)):
                return None
            loops.append(circle)
        outer = max(loops, key=lambda c: c.Radius)
        centre, axis = outer.Center, outer.Axis
        if any(c.Center != centre or c.Axis.cross(axis).Length > 0 for c in loops):
            return None
        wires = [Part.Wire(Part.makeCircle(c.Radius, centre, axis))
                 for c in sorted(loops, key=lambda c: -c.Radius)]
        disk = Part.Face(wires, "Part::FaceMakerBullseye")
        if len(disk.Faces) == 1 and disk.isValid() and disk.Area > 0:
            return disk.Faces[0], centre, axis
    return None


def _phase_bound(face, ideal, origin, axis):
    # Rotation about this line leaves the entire *ideal* point set identical.
    # Choose its seam independently for each trim, then certify the pairing.
    def radial(p):
        return p - origin - axis * (p-origin).dot(axis)
    seam = radial(ideal.Vertexes[0].Point)
    for vertex in face.Vertexes:
        other = radial(vertex.Point)
        if other.Length == 0:
            continue
        angle = math.degrees(math.atan2(axis.dot(seam.cross(other)), seam.dot(other)))
        rotated = ideal.copy()
        rotated.rotate(origin, axis, angle)
        error = _chart_bound(face,rotated)
        if error > TOL:
            error = _rational_bound(face,rotated)
        if error <= TOL:
            return error
    return math.inf


def _radial_rim(wire, centre, axis, radius):
    """Bound every point of a closed profile and prove that it covers all azimuths."""
    if not wire.isClosed():
        return math.inf
    disk = Part.Face(wire)
    if not disk.isValid() or not isinstance(disk.Surface,Part.Plane):
        return math.inf
    plane = disk.Surface
    denom = plane.Axis.dot(axis)
    if not denom:
        return math.inf
    pierce = centre + axis*((plane.Position-centre).dot(plane.Axis)/denom)
    if not disk.isPartOfDomain(*plane.parameter(pierce)):
        return math.inf
    # The axis pierces the interior of a simple planar profile. Every ray
    # in its orthogonal projection therefore meets this complete closed rim.
    rerr = axial = 0.
    for edge in wire.Edges:
        curve = edge.Curve
        if not isinstance(curve,(Part.Circle,Part.Ellipse,Part.BSplineCurve,Part.BezierCurve,Part.Line)):
            return math.inf
        for arc in curve.toBSpline(*edge.ParameterRange).toBezier():
            p = np.array([list(pt-centre) for pt in arc.getPoles()])
            z = p @ np.array(list(axis))
            xy = p - z[:,None]*np.array(list(axis))
            # In the Bernstein product basis the positive denominator and
            # squared radial numerator share the same binomial normalizers.
            w = np.array(arc.getWeights()) * [math.comb(arc.Degree,k) for k in range(len(p))]
            if min(w) <= 0:
                return math.inf
            ww = np.convolve(w,w)
            rr = sum(np.convolve(w*q,w*q) for q in xy.T)/ww
            rerr = max(rerr,abs(math.sqrt(max(0.,min(rr)))-radius),
                       abs(math.sqrt(max(0.,max(rr)))-radius))
            axial = max(axial, max(abs(z)))
    return math.hypot(rerr,axial)


def _prism_bound(template, other):
    """Bound parallel sections throughout a trimmed face of straight extrusion."""
    nominal = template.Surface
    if not isinstance(nominal,(Part.Cylinder,Part.Plane)):
        return math.inf
    chart = _round_chart(template)
    if chart is None:
        return math.inf
    ideal,centre,axis = chart
    if isinstance(nominal,Part.Plane):
        if not isinstance(other.Surface,Part.Plane) or not other.isValid() or len(other.Wires)!=len(ideal.Wires):
            return math.inf
        def circles(face):
            return sorted(face.Wires,key=lambda w:not w.isSame(face.OuterWire))
        area_bound = [_radial_rim(w,centre,axis,rim.Edges[0].Curve.Radius)
                      for w,rim in zip(circles(other),circles(ideal))]
        # The projected trims separate the inner and outer radius bands.
        radii = sorted(e.Curve.Radius for e in ideal.Edges)
        if len(radii) == 2 and 2*max(area_bound) >= radii[1]-radii[0]:
            return math.inf
        # Height and radial extremes may occur at different points on a
        # planar trim. Use their sum to bound the interior as well.
        return sum(area_bound)
    surf = other.Surface
    if not isinstance(surf,Part.SurfaceOfExtrusion) or not isinstance(surf.BasisCurve,Part.BSplineCurve):
        return math.inf
    ref_trim = _phase_bound(template,ideal,centre,axis)
    if ref_trim > TOL:
        return math.inf
    u,U,v,V = other.ParameterRange
    whole = surf.toShape(u,U,v,V)
    gap = _trim_bound(other,whole,match=True)
    if gap is None:
        return math.inf
    spline = surf.BasisCurve
    poles = np.array([list(p) for p in spline.getPoles()])
    weights = np.array(spline.getWeights())
    if min(weights) <= 0:
        return math.inf
    # The derivative of a homogeneous B-spline uses degree / knot span.
    # Translating the numerator to the first pole limits both rational terms.
    d = spline.Degree
    span = min(np.diff(spline.getKnots()))
    diam = max(np.linalg.norm(poles-poles[0],axis=1))
    speed = 3*d*diam*max(weights)/(span*min(weights))
    trim_error = speed*gap[0] + gap[1]
    if trim_error+ref_trim > TOL:
        return math.inf
    rim = Part.Wire(surf.vIso(v).toShape(u,U))
    step = surf.value(u,V)-surf.value(u,v)
    ends = template.ParameterRange[2:]
    dz = ends[1]-ends[0]
    def compare(c,n):
        cap_gap = _radial_rim(rim,c,n,nominal.Radius)
        return cap_gap+(step-n*dz).Length+trim_error+ref_trim
    # From either profile, match an azimuth on the other first section and
    # follow both straight generatrices by the same fraction of their length.
    return min(compare(centre,axis),compare(centre+axis*dz,-axis))


def _trimmed_point_gap(vertex, target, tolerance):
    """An OCCT extremum is an upper estimate; disjointness decides rejection."""
    d = vertex.distToShape(target)[0]
    if d <= tolerance:
        return d
    p = vertex.Point
    for face in target.Faces:
        try:
            uv = face.Surface.parameter(p)
            if face.isPartOfDomain(*uv):
                d = min(d, (face.valueAt(*uv)-p).Length)
                if d <= tolerance:
                    return d
        except (ValueError, RuntimeError, Part.OCCError):
            continue
    # A nonempty common is itself a witness within the ball. A false extremum
    # cannot turn into a rejection unless the trimmed skin misses the ball.
    common = target.common(Part.makeSphere(tolerance,p))
    if common.isNull() or not (common.Faces or common.Edges or common.Vertexes):
        # The disjointness result is a threshold verdict, not a measured
        # distance from the kernel extremum.
        return math.nextafter(tolerance,math.inf)
    return min([tolerance,*[(v.Point-p).Length for v in common.Vertexes]])


def _round_bound(f,g):
    if isinstance(f.Surface, Part.Plane) and isinstance(g.Surface, Part.Plane):
        first, second = _round_chart(f), _round_chart(g)
        if first is not None and second is not None and len(f.Wires) == len(g.Wires):
            radii = [sorted(e.Curve.Radius for e in shape.Edges) for shape,_,_ in (first,second)]
            _, ca, ax = first
            _, cb, bx = second
            # These are complete coaxial circular loops on valid planar faces.
            # Map each radius interval linearly and rotate the plane normals.
            turn = min((ax-bx).Length,(ax+bx).Length)
            return ((ca-cb).Length + max(abs(i-j) for i,j in zip(*radii))
                    + max(map(max,radii))*turn)
    if isinstance(f.Surface, Part.SurfaceOfExtrusion) or isinstance(g.Surface, Part.SurfaceOfExtrusion) or (
            isinstance(f.Surface, Part.Plane) and isinstance(g.Surface, Part.Plane)):
        for template,other in ((f,g),(g,f)):
            error = _prism_bound(template,other)
            if error <= TOL:
                return error
    template, other = f, g
    chart = _round_chart(template)
    if chart is None:
        template, other = g, f
        chart = _round_chart(template)
    if chart is None:
        return math.inf
    ideal, origin, axis = chart
    left = _phase_bound(template, ideal, origin, axis)
    if left > TOL:
        left = _phase_bound(template, ideal.reversed(), origin, axis)
    if left > TOL:
        return math.inf
    # Any rotation Q about the ideal's axis has d(template,Q template)
    # <= d(template,ideal) + d(Q ideal,Q template) = 2*left.
    return 2*left + _phase_bound(other, template, origin, axis)


def _adapt(face):
    """Access OCCT's analytic patch bounds without tessellation."""
    from io import BytesIO
    from OCP.BRep import BRep_Builder
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.BRepTools import BRepTools
    from OCP.TopoDS import TopoDS, TopoDS_Shape
    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, BytesIO(face.exportBrepToString().encode()), BRep_Builder())
    return BRepAdaptor_Surface(TopoDS.Face(shape))


def _bezier_interval(net, first, last, axis):
    """Homogeneous Bernstein net on a subinterval, including endpoint extension."""
    row = np.swapaxes(net,0,axis)
    def split(values,t):
        left,right = [values[0]],[values[-1]]
        while len(values)>1:
            values = (1-t)*values[:-1]+t*values[1:]
            left.append(values[0])
            right.append(values[-1])
        return np.array(left),np.array(right[::-1])
    if first == 1:
        head,_ = split(row,last)
        _,interval = split(head, first/last)
    else:
        _,tail = split(row,first)
        interval,_ = split(tail, (last-first)/(1-first))
    return np.swapaxes(interval,0,axis)


def _bezier_nets(surface, uv):
    """Clip the active knot spans using original (possibly extrapolated) UVs."""
    if isinstance(surface,Part.BezierSurface):
        surface = surface.toBSpline()
    def spans(lo,hi,knots,periodic):
        if periodic and (lo < knots[0] or hi > knots[-1]):
            period = knots[-1]-knots[0]
            start = knots[0] + ((lo-knots[0]) % period)
            end = min(start+hi-lo, start+period)
            pieces = [(start,min(end,knots[-1]))]
            if end > knots[-1]:
                pieces.append((knots[0],end-period))
        else:
            pieces = [(lo,hi)]
        for x,X in pieces:
            for i,(a,b) in enumerate(zip(knots,knots[1:])):
                first = x if i == 0 else max(a,x)
                last = X if i == len(knots)-2 else min(b,X)
                if first < last:
                    yield a,b,(first-a)/(b-a),(last-a)/(b-a)
    u,U,v,V = uv
    us = spans(u,U,surface.getUKnots(),surface.isUPeriodic())
    vs = list(spans(v,V,surface.getVKnots(),surface.isVPeriodic()))
    for a,b,ut,uT in us:
        for c,d,vt,vT in vs:
            span = surface.copy()
            span.segment(a,b,c,d)
            p,w = _net(span)
            net = np.concatenate((p*w[:,:,None], w[:,:,None]),axis=2)
            net = _bezier_interval(_bezier_interval(net,ut,uT,0),vt,vT,1)
            weights = net[:,:,3]
            if weights.min() <= 0:
                raise ValueError("A rational surface patch has no positive denominator bound")
            yield net[:,:,:3]/weights[:,:,None], weights


def _spline_rect(surface, uv):
    return [App.Vector(*point) for p,_ in _bezier_nets(surface,uv) for point in p.reshape(-1,3)]


def _corners(lo, hi):
    return [App.Vector(*p) for p in product(*zip(lo, hi))]


def _offset_spline(surface, uv):
    """Enclose B + offset * normalize(Bu × Bv) over the original UV patch."""
    return [corner for p,w in _bezier_nets(surface.BasisSurface,uv)
            for corner in _offset_net(p,w,surface.OffsetValue)]


def _offset_net(p,w,offset):
    bounds = np.min(p, axis=(0,1)), np.max(p, axis=(0,1))
    # Translate homogeneous numerators to the local hull to avoid dependency
    # on large world coordinates. W is positive; W times each partial
    # derivative equals H' - (B-origin) W' and has the same direction.
    coords = p - bounds[0]
    h = coords*w[:,:,None]
    tangents = []
    for axis in (0,1):
        degree = w.shape[axis]-1
        dh = degree*np.diff(h,axis=axis)
        dw = degree*np.diff(w,axis=axis)
        span = bounds[1]-bounds[0]
        lo = dh.min(axis=(0,1)) - span*max(0.,dw.max())
        hi = dh.max(axis=(0,1)) - span*min(0.,dw.min())
        tangents.append(np.array(list(product(*zip(lo,hi)))))
    # Each component of a cross product is bilinear, with extrema over
    # two boxes at pairs of corners.
    crosses = np.cross(tangents[0][:,None,:], tangents[1][None,:,:])
    lo, hi = crosses.min(axis=(0,1)), crosses.max(axis=(0,1))
    near = np.linalg.norm(np.maximum(0.,np.maximum(lo,-hi)))
    far = np.linalg.norm(np.maximum(np.abs(lo),np.abs(hi)))
    if near > 0:
        normal_lo = np.maximum(-1.,np.minimum(lo/near,lo/far))
        normal_hi = np.minimum(1.,np.maximum(hi/near,hi/far))
    else:
        normal_lo, normal_hi = -np.ones(3), np.ones(3)
    a,b = normal_lo*offset, normal_hi*offset
    return _corners(bounds[0]+np.minimum(a,b), bounds[1]+np.maximum(a,b))


def _canonical_offset(face):
    """Return OCCT's exact analytic offset support, in the original UV chart."""
    from OCP.BRep import BRep_Tool
    from OCP.Geom import Geom_OffsetSurface, Geom_RectangularTrimmedSurface
    surface = BRep_Tool.Surface_s(_adapt(face).Face())
    while isinstance(surface, Geom_RectangularTrimmedSurface):
        surface = surface.BasisSurface()
    return surface.Surface() if isinstance(surface, Geom_OffsetSurface) else None


def _rational_support(surface):
    """These elementary and swept surfaces have exact rational representations."""
    if isinstance(surface, (Part.Plane, Part.Cylinder, Part.Cone, Part.Sphere,
                            Part.Toroid, Part.BezierSurface, Part.BSplineSurface)):
        return True
    if isinstance(surface, Part.RectangularTrimmedSurface):
        return _rational_support(_untrim(surface))
    return isinstance(surface, (Part.SurfaceOfExtrusion, Part.SurfaceOfRevolution)) and isinstance(
        surface.BasisCurve, (Part.Line, Part.Circle, Part.Ellipse, Part.Parabola,
                             Part.Hyperbola, Part.BezierCurve, Part.BSplineCurve))


def _offset_patch(surface, uv):
    # Changing parameterization does not change the point set of this patch.
    # Enclose its rational basis patch and its unit normals on the same patch.
    # Bound both sides if the conversion reverses its parameter orientation.
    patch = surface.BasisSurface.toShape(*uv).toNurbs().Faces[0]
    spline = patch.Surface
    um = sum(uv[:2])/2
    vm = sum(uv[2:])/2
    p = surface.BasisSurface.value(um, vm)
    u, v = spline.parameter(p)
    direction = spline.normal(u, v).dot(surface.BasisSurface.normal(um, vm))
    offset = surface.OffsetValue if direction >= 0 else -surface.OffsetValue
    return [point for poles, weights in _bezier_nets(spline, patch.ParameterRange)
            for point in _offset_net(poles, weights, offset)]


def _enclosure(face):
    from OCP.Bnd import Bnd_Box
    from OCP.BndLib import BndLib_AddSurface
    surface = face.Surface
    box_for_patch = BndLib_AddSurface.Add_s
    if isinstance(surface, (Part.BezierSurface, Part.BSplineSurface)):
        return lambda uv: _spline_rect(surface, uv)
    if isinstance(surface, Part.OffsetSurface):
        if isinstance(surface.BasisSurface, (Part.BezierSurface,Part.BSplineSurface)):
            return lambda uv: _offset_spline(surface, uv)
        # OCCT bounds an offset via its entire offset distance around the
        # basis box. That enclosure does not converge. Its exact elementary
        # equivalent instead admits arbitrarily small boxes in the same UVs.
        from OCP.GeomAdaptor import GeomAdaptor_Surface
        exact = _canonical_offset(face)
        if exact is not None:
            adaptor = GeomAdaptor_Surface(exact)
        elif _rational_support(surface.BasisSurface):
            return lambda uv: _offset_patch(surface, uv)
        else:
            # Ask OCCT for the tight box of the offset patch itself. The
            # ordinary box expands the entire basis by the offset magnitude.
            adaptor = _adapt(face)
            box_for_patch = BndLib_AddSurface.AddOptimal_s
    else:
        adaptor = _adapt(face)
    if isinstance(surface, (Part.SurfaceOfExtrusion,Part.SurfaceOfRevolution)):
        # GeomConvert gives an exact rational swept surface for these generators.
        if isinstance(surface.BasisCurve, (Part.Line,Part.Circle,Part.Ellipse,
                                          Part.Parabola,Part.Hyperbola,
                                          Part.BezierCurve,Part.BSplineCurve)):
            return None, face.toNurbs().Faces[0]
    # Use OCCT's enclosure of the original support. Numerical conversion of
    # a general offset to a spline would change the measured geometry.
    def bound(uv):
        box = Bnd_Box()
        box_for_patch(adaptor, *uv, 0., box)
        coords = box.Get()
        return _corners(coords[:3], coords[3:])
    return bound


def _convex_planar_regions(face):
    """A straight-sided convex planar trim is its own convex hull."""
    if not isinstance(_elementary_basis(face), Part.Plane) or len(face.Wires) != 1:
        return False
    edges = face.OuterWire.OrderedEdges
    if not edges or not all(isinstance(e.Curve, Part.Line) for e in edges):
        return False
    points = [v.Point for v in face.OuterWire.OrderedVertexes]
    n = _elementary_basis(face).Axis
    turns = [(points[i]-points[i-1]).cross(points[(i+1) % len(points)]-points[i]).dot(n)
             for i in range(len(points))]
    return min(turns) >= 0 or max(turns) <= 0


def _planar_cover(face):
    basis = _elementary_basis(face)
    if not isinstance(basis, Part.Plane):
        return None
    return Part.Plane(face.valueAt(0,0), basis.Axis), Part.makeCompound(face.Edges), face


def _inside_planar_trim(cover, poles, radius, mid, tolerance):
    # If the projected control hull fits in a boundary-free disk around an
    # interior projection, the entire patch projects into this trimmed face.
    plane, edges, face = cover
    normal, origin = plane.Axis, plane.Position
    height = max(abs((p-origin).dot(normal)) for p in poles)
    if height > tolerance:
        return math.inf
    centre = mid - normal * (mid-origin).dot(normal)
    if not face.isPartOfDomain(*face.Surface.parameter(centre)):
        return math.inf
    if Part.Vertex(centre).distToShape(edges)[0] >= radius:
        return height
    return math.inf


def _trim_rectangles(face):
    """Enclose each complete trimming arc in its UV control polygon's box."""
    boxes = []
    for wire in face.Wires:
        for edge in wire.OrderedEdges:
            try:
                arcs = _trim_arc(face,edge)
            except (TypeError, ValueError, RuntimeError, Part.OCCError):
                return None
            if not arcs:
                return None
            for poles, weights in arcs:
                if min(weights) <= 0:
                    return None
                lo, hi = poles.min(axis=0), poles.max(axis=0)
                boxes.append((lo[0],hi[0],lo[1],hi[1]))
    return boxes


def _cylinder_frame(face, v):
    x = face.valueAt(0, v)
    c = (x + face.valueAt(math.pi, v)) * .5
    return c, x-c, face.valueAt(math.pi/2,v)-c


def _cylinder_map(source, target):
    """An isometry of angular and axial coordinates, with a global error bound."""
    c,x,y = _cylinder_frame(source,0)
    d,p,q = _cylinder_frame(target,0)
    a = _cylinder_frame(source,1)[0]-c
    b = _cylinder_frame(target,1)[0]-d
    su = 1 if x.cross(y).dot(p.cross(q)) >= 0 else -1
    sv = 1 if a.dot(b) >= 0 else -1
    phase = math.atan2(x.dot(q), x.dot(p))
    offset = (c-d).dot(b) / b.dot(b)
    error = 0.
    for v in source.ParameterRange[2:]:
        t,tx,ty = _cylinder_frame(target,sv*v+offset)
        aligned = (t, math.cos(phase)*tx + math.sin(phase)*ty,
                   su*(-math.sin(phase)*tx + math.cos(phase)*ty))
        # The coefficient differences are affine in v. The triangle inequality
        # bounds the distance on the complete periodic chart between the ends.
        error = max(error, sum((i-j).Length for i,j in zip(_cylinder_frame(source,v), aligned)))
    return su, sv, phase, offset, error


def _cylinder_cover(face, target):
    """Certify a cylindrical face against the union of compatible trimmed charts."""
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.gp import gp_Pnt2d
    from OCP.TopAbs import TopAbs_IN, TopAbs_OUT
    base = face.ParameterRange
    original = [tuple(base)]
    if not isinstance(_elementary_basis(face), Part.Cylinder):
        return 0., original
    src = _trim_rectangles(face)
    if src is None:
        return 0., original

    # Transform complete pcurve boxes into the source chart. A connected open
    # cell untouched by trims has one winding classification across its interior.
    period = 2*math.pi
    copies = []
    boundaries = list(src)
    carrier_gap = 0.
    for other in target.Faces:
        if not isinstance(_elementary_basis(other), Part.Cylinder):
            continue
        su, sv, phase, offset, error = _cylinder_map(face,other)
        if error > TOL:
            continue
        dst = _trim_rectangles(other)
        if dst is None:
            continue
        def in_source(box):
            u,U,v,V = box
            return (*sorted(((u-phase)*su,(U-phase)*su)),
                    *sorted(((v-offset)*sv,(V-offset)*sv)))
        t0,t1,w0,w1 = in_source(other.ParameterRange)
        k0 = math.ceil((base[0]-t1)/period)
        k1 = math.floor((base[1]-t0)/period)
        tf = _adapt(other).Face()
        for k in range(k0,k1+1):
            shift = k*period
            domain = (t0+shift,t1+shift,w0,w1)
            boxes = [(u+shift,U+shift,v,V) for u,U,v,V in map(in_source,dst)]
            boundaries.extend([domain,*boxes])
            copies.append((domain,boxes,tf,(su,sv,phase-su*shift,offset)))
            carrier_gap = max(carrier_gap,error)
    if not copies:
        return 0.,original
    def cuts(start,end,axis):
        return sorted({start,end} | {x for b in boundaries for x in b[axis:axis+2]
                                     if start < x < end})
    us,vs = cuts(*base[:2],0),cuts(*base[2:],2)
    if (len(us)-1)*(len(vs)-1) > 256:
        # A dense trim arrangement is cheaper to handle by spatial refinement.
        return 0., original
    sf = _adapt(face).Face()
    def state(cell, boxes, shape, change):
        u,U,v,V = cell
        if any(u<h and l<U and v<H and L<V for l,h,L,H in boxes):
            return None
        su,sv,phase,offset = change
        mid = gp_Pnt2d(su*(u+U)/2+phase, sv*(v+V)/2+offset)
        return BRepClass_FaceClassifier(shape,mid,0.).State()

    certified, pending = [], []
    for u,U in zip(us,us[1:]):
        for v,V in zip(vs,vs[1:]):
            cell = u,U,v,V
            if any(lo<=u and U<=hi and low<=v and V<=high
                   and state(cell,boxes,tf,change) == TopAbs_IN
                   for (lo,hi,low,high),boxes,tf,change in copies):
                certified.append(cell)
            elif state(cell,src,sf,(1,1,0.,0.)) != TopAbs_OUT:
                pending.append(cell)
    if not certified:
        return 0., original

    # Include slender cells by clamping each parameter into a certified
    # neighboring trim. Cylinder surface speed is radius in u and 1 in v.
    radius = _cylinder_frame(face,0)[1].Length
    residual = []
    upper = carrier_gap
    for u,U,v,V in pending:
        d = min(radius*max(0.,l-u,U-h) + max(0.,L-v,V-H)
                for l,h,L,H in certified) + carrier_gap
        if d <= TOL:
            upper = max(upper,d)
        else:
            residual.append((u,U,v,V))
    return upper,residual


def _contained_in_normal_prisms(face, target, tolerance):
    """Try a whole-face BRep inclusion proof using subsets of the target's tolerance skin."""
    from OCP.Precision import Precision

    # Each point in a prism is within its half-height of its actual trimmed
    # base. The round tools likewise have a point on an edge or vertex at
    # most 'radius' away. Leave room for the stored OCCT Boolean tolerances.
    reserve = max(Precision.Confusion_s(), face.getTolerance(1), target.getTolerance(1))
    radius = tolerance - 8*reserve
    if radius <= 0:
        return False
    search = face.BoundBox
    search.enlarge(tolerance)
    def tools():
        for patch in target.Faces:
            if isinstance(patch.Surface, Part.Plane) and search.intersect(patch.BoundBox):
                dz = radius * patch.Surface.Axis
                lo = patch.copy()
                lo.translate(-dz)
                yield lo.extrude(2*dz)
        for edge in target.Edges:
            if not search.intersect(edge.BoundBox):
                continue
            curve = edge.Curve
            if isinstance(curve, Part.Line) and len(edge.Vertexes) == 2:
                start, end = (v.Point for v in edge.Vertexes)
                yield Part.makeCylinder(radius, (end-start).Length, start, end-start)
            elif (isinstance(curve, Part.Circle) and edge.isClosed()
                  and math.isclose(edge.LastParameter-edge.FirstParameter, math.tau,
                                   rel_tol=32*math.ulp(1.)) and curve.Radius > radius):
                yield Part.makeTorus(curve.Radius, radius, curve.Center, curve.Axis)
        for vertex in target.Vertexes:
            if search.intersect(vertex.BoundBox):
                yield Part.makeSphere(radius, vertex.Point)

    remaining = face
    for solid in tools():
        if not remaining.BoundBox.intersect(solid.BoundBox) or not solid.isValid():
            continue
        if solid.getTolerance(1) > 2*reserve:
            continue
        rest = remaining.cut(solid)
        if rest.isNull() or not rest.isValid():
            return False
        if not rest.Vertexes:
            return True
        if rest.getTolerance(1) > 2*reserve:
            return False
        remaining = rest
    return False


def _adaptive(face, target, convex, tolerance, work):
    """Cover a trimmed face with certified spatial enclosures until decided."""
    enclose = _enclosure(face)
    if isinstance(enclose, tuple):
        return _adaptive(enclose[1], target, convex, tolerance, work)
    planes = [r for f in target.Faces if (r := _planar_cover(f)) is not None]
    surface = face.Surface
    base = face.ParameterRange
    us, vs = [base[0],base[1]], [base[2],base[3]]
    if isinstance(surface, Part.BSplineSurface):
        # Knot spans expose highly localized detail at the first refinement.
        us = sorted(set(us + [u for u in surface.getUKnots() if us[0] < u < us[-1]]))
        vs = sorted(set(vs + [v for v in surface.getVKnots() if vs[0] < v < vs[-1]]))
    pending = [(u,U,v,V) for u,U in zip(us, us[1:]) for v,V in zip(vs, vs[1:])]
    if isinstance(_elementary_basis(face), Part.Cylinder):
        max_bound, pending = _cylinder_cover(face,target)
    else:
        max_bound = 0.
    while pending:
        if next(work, None) is None:
            if _contained_in_normal_prisms(face, target, tolerance):
                return tolerance
            raise ArithmeticError("Boundary comparison inconclusive: trimmed charts, spatial "
                                  "enclosures and BRep containment did not certify this face")
        uv = pending.pop()
        u,U,v,V = uv
        um,vm = (u+U)*.5,(v+V)*.5
        mid = face.valueAt(um,vm)
        if face.isPartOfDomain(um,vm):
            d = _trimmed_point_gap(Part.Vertex(mid), target, tolerance)
            if d > tolerance:
                return d
        poles = enclose(uv)
        # The complete patch lies in the convex hull of its poles/box corners.
        # Distance to any convex subset of the other skin is convex.
        for region in convex:
            upper = max(Part.Vertex(p).distToShape(region)[0] for p in poles)
            if upper <= tolerance:
                max_bound = max(max_bound, upper)
                break
        else:
            radius = max((p-mid).Length for p in poles)
            planar = min((_inside_planar_trim(cover, poles, radius, mid, tolerance)
                          for cover in planes), default=math.inf)
            if planar <= tolerance:
                max_bound = max(max_bound, planar)
                continue
            point = Part.Vertex(mid)
            if not face.isPartOfDomain(um,vm) and radius and _trimmed_point_gap(point,face,radius) > radius:
                continue
            d = _trimmed_point_gap(point,target,tolerance)
            if d+radius <= tolerance:
                max_bound = max(max_bound, d+radius)
                continue
            # A full periodic span can have identical endpoints. Its two
            # half chords still expose the direction that needs subdivision.
            ul = (face.valueAt(U, vm)-mid).Length + (face.valueAt(u, vm)-mid).Length
            vl = (face.valueAt(um, V)-mid).Length + (face.valueAt(um, v)-mid).Length
            split_u, split_v = u < um < U, v < vm < V
            if split_u and (not split_v or ul >= vl):
                pending.extend([(u,um,v,V), (um,U,v,V)])
            elif split_v:
                pending.extend([(u,U,v,vm), (u,U,vm,V)])
            else:
                # At adjacent floating-point UVs a generic kernel box may
                # still contain a whole control hull. Recompute the enclosure
                # on the clipped rational chart, whose range converges.
                if _rational_support(surface):
                    chart = surface.toShape(*uv).toNurbs().Faces[0]
                    support = _spline_rect(chart.Surface, chart.ParameterRange)
                    radius = max((p-mid).Length for p in support)
                if d+radius > tolerance:
                    raise ArithmeticError("Could not resolve the boundary-distance interval")
                max_bound = max(max_bound, d+radius)
    return max_bound


def _planar_unions(shape):
    """Merge face partitions only when their exposed trimming edges survive."""
    planar = [f for f in shape.Faces if isinstance(f.Surface, Part.Plane)]
    pieces = [f for f in shape.Faces if not isinstance(f.Surface, Part.Plane)]
    error, changed = 0., False
    while planar:
        first = planar.pop()
        plane = first.Surface
        group = [first]
        for f in planar[:]:
            if f.Surface.Axis.cross(plane.Axis).Length * shape.BoundBox.DiagonalLength > TOL:
                continue
            if any(abs((v.Point-plane.Position).dot(plane.Axis)) > TOL for v in f.Vertexes):
                continue
            group.append(f)
            planar.remove(f)
        if len(group) == 1:
            pieces.extend(group)
            continue
        merged = Part.makeShell(group).removeSplitter().Faces
        if len(merged) != 1 or not merged[0].isValid() or not isinstance(merged[0].Surface, Part.Plane):
            pieces.extend(group)
            continue
        # Shared edges cancel in the planar projection. Equality of the
        # remaining physical trimming edges establishes identical projected
        # domains (parity), including the holes of non-convex faces.
        all_edges = [e for f in group for e in f.Edges]
        rim = [e for e in all_edges if sum(e.isSame(k) for k in all_edges) == 1]
        if (len(rim) != len(merged[0].Edges)
                or not all(any(e.isSame(k) for k in merged[0].Edges) for e in rim)):
            pieces.extend(group)
            continue
        target = merged[0].Surface
        # Bound support-plane separation on each complete trimmed face.
        delta = max(abs((f.valueAt(u,v)-target.Position).dot(target.Axis))
                    for f in group for u,v in product(f.ParameterRange[:2],f.ParameterRange[2:]))
        pieces.extend(merged)
        changed = True
        error = max(error, delta)
    return (Part.makeCompound(pieces), error) if changed else (shape, 0.)


def boundary_error(a, b):
    """Certify complete boundaries in both directions, including trimmed edges."""
    rigid = _rigid(a, b)
    if rigid <= TOL:
        return rigid
    # A bound on refinement work, independent of the program being graded.
    # Both directions and all faces of this pair draw from the same budget.
    work = iter(range(4096))
    if len(a.Faces) != len(b.Faces):
        first, da = _planar_unions(a)
        second, db = _planar_unions(b)
        if first is not a or second is not b:
            try:
                merged = _boundary_charts(first, second, work)
            except ArithmeticError:
                merged = math.inf
            if merged + da + db <= TOL:
                return merged + da + db
    return _boundary_charts(a, b, work)


def _boundary_charts(a, b, work):
    # Paired analytic charts are only an acceleration. Reordered/split/rebuilt
    # faces can be covered by the entire opposing boundary in either direction.
    covered_a, covered_b = set(), set()
    worst = 0.
    for i, (f, g) in enumerate(zip(a.Faces, b.Faces)):
        error = _chart_bound(f, g)
        if error <= TOL:
            worst = max(worst, error)
            covered_a.add(i)
            covered_b.add(i)
    af, bf = a.Faces, b.Faces
    for i, f in enumerate(af):
        if i in covered_a:
            continue
        options = sorted((j for j in range(len(bf)) if j not in covered_b
                          and type(_elementary_basis(f)) == type(_elementary_basis(bf[j]))),
                         key=lambda j: (f.CenterOfMass-bf[j].CenterOfMass).Length)
        for j in options:
            error = _chart_bound(f,bf[j])
            if error <= TOL:
                worst = max(worst,error)
                covered_a.add(i)
                covered_b.add(j)
                break
    for i, f in enumerate(af):
        if i in covered_a:
            continue
        options = sorted((j for j in range(len(bf)) if j not in covered_b),
                         key=lambda j: (f.CenterOfMass-bf[j].CenterOfMass).Length)
        for j in options:
            try:
                error = _rational_bound(f,bf[j])
                if error > TOL:
                    error = _round_bound(f,bf[j])
            except (ValueError, RuntimeError, Part.OCCError):
                # A degenerate chart can still be covered in physical space.
                continue
            if error <= TOL:
                worst = max(worst,error)
                covered_a.add(i)
                covered_b.add(j)
                break
    remaining = []
    for source, other, done in ((a, b, covered_a), (b, a, covered_b)):
        faces = [f for i,f in enumerate(source.Faces) if i not in done]
        if not faces:
            continue
        skin = Part.makeCompound(other.Faces)
        convex = [f for f in other.Faces if _convex_planar_regions(f)]
        for face in faces:
            # Exact points provide early rejection only. A spline span can be
            # extremely small relative to the parent face's parameter range.
            for p in face.Vertexes:
                d = _trimmed_point_gap(p,skin,TOL)
                if d > TOL:
                    return d
            surf = face.Surface
            if isinstance(surf, Part.BSplineSurface):
                us, vs = surf.getUKnots(),surf.getVKnots()
                for u,U in zip(us,us[1:]):
                    for v,V in zip(vs,vs[1:]):
                        um,vm = (u+U)*.5,(v+V)*.5
                        if face.isPartOfDomain(um,vm):
                            d = _trimmed_point_gap(Part.Vertex(face.valueAt(um,vm)),skin,TOL)
                            if d > TOL:
                                return d
            remaining.append((face,skin,convex))
    for face,skin,convex in remaining:
        error = _adaptive(face, skin, convex, TOL, work)
        if error > TOL:
            return error
        worst = max(worst, error)
    return worst
