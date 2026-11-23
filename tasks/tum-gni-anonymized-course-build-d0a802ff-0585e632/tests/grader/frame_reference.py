"""Plane-frame reference for the carrier chart of the roof trusses.

Everything is in N and mm.  The model follows the brief literally:

* each truss is a vertical plane frame on the placed Axis segments;
* both chords are continuous Euler-Bernoulli beam-columns, the webs are
  pin-ended bars with their analytical length;
* sections are integrated exactly from the swept profile of each placed
  Body (straight edges and circular arcs), keeping the profile roll;
* the lower chord is pinned at its station origin and rests on a
  vertical-reaction roller at the far end;
* self-weight is uniform on the chords and lumped half-and-half at the
  web ends. The wheels and the suspended unit act on the lower chord axis.

The frame is solved by the direct stiffness method with the joints as
nodes.  Loads inside a chord panel enter through their fixed-end actions,
which are obtained from the clamped two-point problem of the element
state (w, w', EI w'', EI w''') and not from tabulated formulas; the same
state is integrated along each panel afterwards, so moments and
deflections between the joints are exact.
"""
import math

import numpy as np
import scipy.linalg

E_STEEL = 200000.0
UNIT_WEIGHT = 0.0000770
WHEELBASE = 1100.0
SUSPENDED_FORCE = 40000.0
CLEARANCE = 400.0
HEAVY_SHARE = 0.65
LIGHT_SHARE = 0.35
STRESS_LIMIT = 160.0
SAG_SPAN_RATIO = 600.0
SAG_CAP = 20.0
X_TIE = 0.1
CLUSTER = 0.5

MODES = ("STRESS", "SAG", "WEB_BUCKLING")


# --------------------------------------------------------------------------
# IFC reading
# --------------------------------------------------------------------------
def _triple(values):
    out = np.zeros(3)
    vals = list(values)
    out[: len(vals)] = vals
    return out


def _placement3d(placement):
    """4x4 matrix of an IfcAxis2Placement3D."""
    mat = np.eye(4)
    if placement is None:
        return mat
    z = _triple(placement.Axis.DirectionRatios) if placement.Axis else np.array([0.0, 0.0, 1.0])
    x = _triple(placement.RefDirection.DirectionRatios) if placement.RefDirection else np.array([1.0, 0.0, 0.0])
    z = z / np.linalg.norm(z)
    x = x - (x @ z) * z
    x = x / np.linalg.norm(x)
    mat[:3, 0] = x
    mat[:3, 1] = np.cross(z, x)
    mat[:3, 2] = z
    mat[:3, 3] = _triple(placement.Location.Coordinates)
    return mat


def _object_placement(local):
    """World matrix of an IfcLocalPlacement chain."""
    mat = np.eye(4)
    while local is not None:
        mat = _placement3d(local.RelativePlacement) @ mat
        local = local.PlacementRelTo
    return mat


def _mapped_matrix(item):
    """Matrix carried by an IfcMappedItem (uniform operators without axes)."""
    target = item.MappingTarget
    if any(getattr(target, name, None) is not None for name in ("Axis1", "Axis2", "Axis3")):
        raise ValueError("mapped item with a rotated operator is not supported")
    if target.is_a("IfcCartesianTransformationOperator3DnonUniform"):
        raise ValueError("non-uniform mapped item is not supported")
    scale = 1.0 if target.Scale is None else float(target.Scale)
    op = np.eye(4)
    op[:3, :3] *= scale
    op[:3, 3] = _triple(target.LocalOrigin.Coordinates)
    origin = _placement3d(item.MappingSource.MappingOrigin)
    return op @ np.linalg.inv(origin)


def _representation(product, identifier):
    for rep in product.Representation.Representations:
        if rep.RepresentationIdentifier == identifier:
            return rep
    raise ValueError("%s has no %s representation" % (product.GlobalId, identifier))


def _leaf_items(rep, base):
    """(item, matrix) pairs of a representation with mapped items resolved."""
    out = []
    for item in rep.Items:
        if item.is_a("IfcMappedItem"):
            out.extend(_leaf_items(item.MappingSource.MappedRepresentation, base @ _mapped_matrix(item)))
        else:
            out.append((item, base))
    return out


def _length_factor(model):
    """Millimetres per model length unit."""
    for unit in model.by_type("IfcUnitAssignment")[0].Units:
        if unit.is_a("IfcSIUnit") and unit.UnitType == "LENGTHUNIT":
            prefix = {None: 1000.0, "MILLI": 1.0, "CENTI": 10.0, "DECI": 100.0}
            return prefix[unit.Prefix]
    raise ValueError("no SI length unit")


def _angle_factor(model):
    """Radians per model plane-angle unit."""
    for unit in model.by_type("IfcUnitAssignment")[0].Units:
        if unit.is_a("IfcConversionBasedUnit") and unit.UnitType == "PLANEANGLEUNIT":
            return float(unit.ConversionFactor.ValueComponent.wrappedValue)
    return 1.0


def _axis_ends(member, world, to_mm):
    """World end points (mm) of the Axis segment of a member."""
    points = []
    for item, mat in _leaf_items(_representation(member, "Axis"), world):
        local = []
        if item.is_a("IfcTrimmedCurve") and item.BasisCurve.is_a("IfcLine"):
            for trim in (item.Trim1, item.Trim2):
                cart = [t for t in trim if hasattr(t, "is_a") and t.is_a("IfcCartesianPoint")]
                if not cart:
                    raise ValueError("axis trims without points")
                local.append(_triple(cart[0].Coordinates))
        elif item.is_a("IfcPolyline"):
            local.extend(_triple(p.Coordinates) for p in item.Points)
        else:
            raise ValueError("unsupported axis item %s" % item.is_a())
        points.extend(mat[:3, :3] @ p + mat[:3, 3] for p in local)
    if len(points) != 2:
        raise ValueError("axis of %s is not a single segment" % member.GlobalId)
    return points[0] * to_mm, points[1] * to_mm


_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(16)


def _profile_pieces(curve, angle):
    """Boundary pieces of a closed profile curve, in traversal order.

    A piece is ("L", start, end) or ("A", centre, radius, angle0, angle1) with
    the arc swept from angle0 to angle1 (radians, signed)."""
    pieces = []
    if curve.is_a("IfcPolyline"):
        pts = [np.array(p.Coordinates[:2], dtype=float) for p in curve.Points]
        for a, b in zip(pts[:-1], pts[1:]):
            pieces.append(("L", a, b))
        if np.linalg.norm(pts[0] - pts[-1]) > 1e-9:
            pieces.append(("L", pts[-1], pts[0]))
        return pieces
    if not curve.is_a("IfcCompositeCurve"):
        raise ValueError("unsupported profile curve %s" % curve.is_a())
    for seg in curve.Segments:
        parent = seg.ParentCurve
        local = []
        if parent.is_a("IfcPolyline"):
            pts = [np.array(p.Coordinates[:2], dtype=float) for p in parent.Points]
            local = [("L", a, b) for a, b in zip(pts[:-1], pts[1:])]
        elif parent.is_a("IfcTrimmedCurve") and parent.BasisCurve.is_a("IfcCircle"):
            circle = parent.BasisCurve
            centre = np.array(circle.Position.Location.Coordinates[:2], dtype=float)
            ref = circle.Position.RefDirection
            phase = math.atan2(ref.DirectionRatios[1], ref.DirectionRatios[0]) if ref else 0.0
            trims = []
            for trim in (parent.Trim1, parent.Trim2):
                par = [t for t in trim if not (hasattr(t, "is_a") and t.is_a("IfcCartesianPoint"))]
                if not par:
                    raise ValueError("arc trims without parameters")
                trims.append(float(par[0].wrappedValue) * angle)
            a0, a1 = trims
            if parent.SenseAgreement:
                while a1 <= a0 + 1e-12:
                    a1 += 2.0 * math.pi
            else:
                while a1 >= a0 - 1e-12:
                    a1 -= 2.0 * math.pi
            local = [("A", centre, float(circle.Radius), a0 + phase, a1 + phase)]
        else:
            raise ValueError("unsupported profile segment %s" % parent.is_a())
        if not seg.SameSense:
            flipped = []
            for piece in reversed(local):
                if piece[0] == "L":
                    flipped.append(("L", piece[2], piece[1]))
                else:
                    flipped.append(("A", piece[1], piece[2], piece[4], piece[3]))
            local = flipped
        pieces.extend(local)
    return pieces


def _piece_point(piece, u):
    if piece[0] == "L":
        return piece[1] + (piece[2] - piece[1]) * u
    ang = piece[3] + (piece[4] - piece[3]) * u
    return piece[1] + piece[2] * np.array([math.cos(ang), math.sin(ang)])


def profile_properties(curve, angle=1.0):
    """Exact area, centroid and central second moments of a closed profile.

    Returns area, centroid (2,), J (2x2 with J[i, j] = integral r_i r_j dA about
    the centroid) and an array of boundary points for extreme-fibre searches."""
    pieces = _profile_pieces(curve, angle)
    for a, b in zip(pieces, pieces[1:] + pieces[:1]):
        if np.linalg.norm(_piece_point(a, 1.0) - _piece_point(b, 0.0)) > 1e-6:
            raise ValueError("profile boundary is not closed")
    area = sx = sy = ixx = iyy = ixy = 0.0
    outline = []
    for piece in pieces:
        if piece[0] == "L":
            spans = [(0.0, 1.0)]
        else:
            count = max(1, int(math.ceil(abs(piece[4] - piece[3]) / (math.pi / 4.0))))
            spans = [(k / count, (k + 1) / count) for k in range(count)]
        for u0, u1 in spans:
            half = 0.5 * (u1 - u0)
            for gx, gw in zip(_GAUSS_X, _GAUSS_W):
                u = u0 + half * (gx + 1.0)
                if piece[0] == "L":
                    x, y = piece[1] + (piece[2] - piece[1]) * u
                    dx, dy = piece[2] - piece[1]
                else:
                    sweep = piece[4] - piece[3]
                    ang = piece[3] + sweep * u
                    x = piece[1][0] + piece[2] * math.cos(ang)
                    y = piece[1][1] + piece[2] * math.sin(ang)
                    dx = -piece[2] * math.sin(ang) * sweep
                    dy = piece[2] * math.cos(ang) * sweep
                w = gw * half
                area += 0.5 * (x * dy - y * dx) * w
                sx += 0.5 * x * x * dy * w
                sy += -0.5 * y * y * dx * w
                iyy += x ** 3 * dy * w / 3.0
                ixx += -(y ** 3) * dx * w / 3.0
                ixy += 0.5 * x * x * y * dy * w
        samples = 2 if piece[0] == "L" else 721
        outline.extend(_piece_point(piece, k / (samples - 1.0)) for k in range(samples))
    if area < 0.0:
        area, sx, sy, ixx, iyy, ixy = -area, -sx, -sy, -ixx, -iyy, -ixy
    centroid = np.array([sx / area, sy / area])
    jxx = iyy - area * centroid[0] ** 2
    jyy = ixx - area * centroid[1] ** 2
    jxy = ixy - area * centroid[0] * centroid[1]
    return area, centroid, np.array([[jxx, jxy], [jxy, jyy]]), np.array(outline)


def _body_section(member, world, mid, axis_dir, to_mm, angle, cache):
    """Section of the placed Body normal to the member axis at mid-length.

    Returns area, the 3x3 tensor of central second moments in world axes, the
    world outline points and the world centroid (all in mm)."""
    leaves = _leaf_items(_representation(member, "Body"), world)
    if len(leaves) != 1 or not leaves[0][0].is_a("IfcExtrudedAreaSolid"):
        raise ValueError("Body of %s is not one extruded solid" % member.GlobalId)
    solid, mat = leaves[0]
    frame = mat @ _placement3d(solid.Position)
    direction = _triple(solid.ExtrudedDirection.DirectionRatios)
    direction = direction / np.linalg.norm(direction)
    if abs(abs(direction[2]) - 1.0) > 1e-9:
        raise ValueError("oblique extrusion is not supported")
    scale = np.linalg.norm(frame[:3, 0])
    sweep = frame[:3, :3] @ direction / scale
    if abs(abs(sweep @ axis_dir) - 1.0) > 1e-9:
        raise ValueError("Body of %s is not swept along its axis" % member.GlobalId)
    start = frame[:3, 3] * to_mm
    reach = (mid - start) @ sweep
    depth = float(solid.Depth) * scale * to_mm
    if not (-1e-6 <= reach <= depth + 1e-6):
        raise ValueError("Body of %s does not reach its mid-length" % member.GlobalId)
    profile = solid.SweptArea
    if not profile.is_a("IfcArbitraryClosedProfileDef") or profile.is_a("IfcArbitraryProfileDefWithVoids"):
        raise ValueError("unsupported profile %s" % profile.is_a())
    if profile.id() not in cache:
        cache[profile.id()] = profile_properties(profile.OuterCurve, angle)
    area, centroid, j2, outline = cache[profile.id()]
    basis = frame[:3, :2] * to_mm          # world image of the profile axes, in mm per profile unit
    area = area * (scale * to_mm) ** 2
    tensor = basis @ j2 @ basis.T * (scale * to_mm) ** 2
    pts = (outline - centroid) @ basis.T
    centre = start + basis @ centroid + sweep * reach
    return area, tensor, pts, centre


# --------------------------------------------------------------------------
# beam element with loads inside the span
# --------------------------------------------------------------------------
def _beam_stiffness(ei, length):
    l = length
    return ei / l ** 3 * np.array([
        [12.0, 6.0 * l, -12.0, 6.0 * l],
        [6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l],
        [-12.0, -6.0 * l, 12.0, -6.0 * l],
        [6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l],
    ])


def _clamped_start(length, w_end, t_end):
    """Start moment and shear (EI w'', EI w''') of a clamped-clamped span.

    w_end and t_end are EI times the end deflection and rotation that the span
    loads alone would produce from a zero start state; the start actions are the
    ones that cancel them."""
    l = length
    m0 = (6.0 * w_end - 2.0 * t_end * l) / l ** 2
    s0 = (6.0 * t_end * l - 12.0 * w_end) / l ** 3
    return m0, s0


def _uniform_fixed_end(length, q):
    """Reactions on a clamped span for a downward uniform load q.

    Returns the four end actions (F_i, M_i, F_j, M_j) acting on the element in
    the sense of its (w, w') degrees of freedom."""
    l = length
    m0, s0 = _clamped_start(l, q * l ** 4 / 24.0, q * l ** 3 / 6.0)
    m_end = m0 + s0 * l - q * l * l / 2.0
    s_end = s0 - q * l
    return np.array([s0, -m0, -s_end, m_end])


class Truss(object):
    """One truss assembly: geometry, sections, frame model and rating."""

    def __init__(self, guid, members):
        self.guid = guid
        self.members = members
        self._build_geometry()
        self._build_frame()

    # ------------------------------------------------------------ geometry
    def _build_geometry(self):
        members = self.members
        flat = [m for m in members if abs(m["p1"][2] - m["p2"][2]) <= 0.05 + 1e-4 * m["length"]]
        if len(flat) != 2:
            raise ValueError("truss %s: expected two horizontal chords, found %d" % (self.guid, len(flat)))
        low, top = sorted(flat, key=lambda m: m["p1"][2] + m["p2"][2])
        p, q = low["p1"], low["p2"]
        if abs(p[0] - q[0]) <= X_TIE:
            origin, far = (p, q) if p[1] <= q[1] else (q, p)
        else:
            origin, far = (p, q) if p[0] < q[0] else (q, p)
        self.low, self.top = low, top
        self.origin, self.far = origin, far
        self.span = float(np.linalg.norm(far - origin))
        ex = far - origin
        ex = ex / np.linalg.norm(ex)
        self.ex = ex
        up = np.array([0.0, 0.0, 1.0])
        self.normal = np.cross(ex, up)
        self.normal /= np.linalg.norm(self.normal)
        z_low = 0.5 * (low["p1"][2] + low["p2"][2])
        z_top = 0.5 * (top["p1"][2] + top["p2"][2])
        self.z_low = z_low
        self.depth = float(z_top - z_low)
        self.s_max = self.span - WHEELBASE
        self.sag_limit = min(self.span / SAG_SPAN_RATIO, SAG_CAP)

        def station(point):
            return float((point - origin) @ ex)

        self.station = station
        t1, t2 = sorted([station(top["p1"]), station(top["p2"])])
        for chord in (low, top):
            chord["kind"] = "chord"
            vertical = np.cross(self.normal, ex)
            vertical /= np.linalg.norm(vertical)
            chord["I"] = float(vertical @ chord["tensor"] @ vertical)
            chord["c"] = float(np.max(np.abs(chord["outline"] @ vertical)))
        webs = [m for m in members if m is not low and m is not top]
        low_st, top_st = [], []
        for web in webs:
            web["kind"] = "web"
            axis = (web["p2"] - web["p1"]) / web["length"]
            e1 = self.normal
            e2 = np.cross(axis, e1)
            e2 /= np.linalg.norm(e2)
            frame = np.column_stack([e1, e2])
            web["I"] = float(np.linalg.eigvalsh(frame.T @ web["tensor"] @ frame)[0])
            web["c"] = None
            ends = []
            for point in (web["p1"], web["p2"]):
                height = point[2] - z_low
                on_low = abs(height) <= abs(height - self.depth)
                ends.append((on_low, station(point)))
                (low_st if on_low else top_st).append(station(point))
            web["ends_plane"] = ends
        self.webs = webs
        self.xb = self._cluster(low_st, 0.0, self.span)
        self.xt = self._cluster(top_st, t1, t2)
        for web in webs:
            nodes = []
            for on_low, st in web["ends_plane"]:
                if on_low:
                    nodes.append(int(np.argmin(np.abs(self.xb - st))))
                else:
                    nodes.append(len(self.xb) + int(np.argmin(np.abs(self.xt - st))))
            web["nodes"] = nodes

    @staticmethod
    def _cluster(stations, first, last):
        values = sorted(stations)
        groups = []
        for value in values:
            if groups and value - groups[-1][-1] <= CLUSTER:
                groups[-1].append(value)
            else:
                groups.append([value])
        reps = [float(np.mean(g)) for g in groups]
        if not reps or abs(reps[0] - first) > CLUSTER:
            reps.insert(0, first)
        if abs(reps[-1] - last) > CLUSTER:
            reps.append(last)
        reps[0], reps[-1] = first, last
        return np.array(reps)

    # --------------------------------------------------------------- frame
    def _build_frame(self):
        xb, xt = self.xb, self.xt
        nb, nt = len(xb) - 1, len(xt) - 1
        self.nb, self.nt = nb, nt
        self.lb, self.lt = np.diff(xb), np.diff(xt)
        nn = nb + nt + 2
        nd = 3 * nn
        node_x = np.concatenate([xb, xt])
        node_z = np.concatenate([np.zeros(nb + 1), np.full(nt + 1, self.depth)])
        K = np.zeros((nd, nd))
        F = np.zeros(nd)
        self.chords = []
        for chord, first, lengths in ((self.low, 0, self.lb), (self.top, nb + 1, self.lt)):
            ea, ei = E_STEEL * chord["A"], E_STEEL * chord["I"]
            q = UNIT_WEIGHT * chord["A"]
            kb, fixed = [], []
            for e, l in enumerate(lengths):
                i, j = first + e, first + e + 1
                ax = np.array([3 * i, 3 * j])
                K[np.ix_(ax, ax)] += ea / l * np.array([[1.0, -1.0], [-1.0, 1.0]])
                bd = np.array([3 * i + 1, 3 * i + 2, 3 * j + 1, 3 * j + 2])
                k = _beam_stiffness(ei, l)
                K[np.ix_(bd, bd)] += k
                r = _uniform_fixed_end(l, q)
                F[bd] -= r
                kb.append(k)
                fixed.append(r)
            self.chords.append({"member": chord, "first": first, "lengths": np.asarray(lengths, dtype=float),
                                "ea": ea, "ei": ei, "q": q, "kb": np.array(kb), "fixed": np.array(fixed),
                                "x0": node_x[first:first + len(lengths)].copy()})
        nw = len(self.webs)
        self.web_rows = np.zeros((nw, nd))
        self.web_area = np.array([w["A"] for w in self.webs])
        self.web_crit = np.array([math.pi ** 2 * E_STEEL * w["I"] / (2.0 * w["length"] ** 2) for w in self.webs])
        for k, web in enumerate(self.webs):
            a, b = web["nodes"]
            dx, dz = node_x[b] - node_x[a], node_z[b] - node_z[a]
            h = math.hypot(dx, dz)
            g = np.zeros(nd)
            g[3 * a], g[3 * a + 1] = -dx / h, -dz / h
            g[3 * b], g[3 * b + 1] = dx / h, dz / h
            stiff = E_STEEL * web["A"] / web["length"]
            K += stiff * np.outer(g, g)
            self.web_rows[k] = stiff * g
            weight = UNIT_WEIGHT * web["A"] * web["length"]
            F[3 * a + 1] -= 0.5 * weight
            F[3 * b + 1] -= 0.5 * weight
        fixed_dofs = [0, 1, 3 * nb + 1]
        free = np.setdiff1d(np.arange(nd), fixed_dofs)
        self.nd, self.free = nd, free
        kff = K[np.ix_(free, free)]
        scale = 1.0 / np.sqrt(np.diag(kff))
        self._scale = scale
        self._factor = scipy.linalg.cho_factor(kff * scale[:, None] * scale[None, :])
        self.u_perm = self._solve(F[:, None])[:, 0]
        self.perm = self._chord_states(self.u_perm[:, None], permanent=True)
        self.web_perm = self.web_rows @ self.u_perm

    def _solve(self, rhs):
        """Displacements for load columns rhs (nd x k)."""
        out = np.zeros_like(rhs)
        scaled = rhs[self.free] * self._scale[:, None]
        out[self.free] = scipy.linalg.cho_solve(self._factor, scaled) * self._scale[:, None]
        return out

    def _chord_states(self, u, permanent, extra=None):
        """Panel start states of both chords for displacement columns u (nd x k).

        Returns, per chord, arrays (k, panels) of the axial force and of the start
        state w, w', EI w'', EI w'''.  `extra` adds fixed-end actions of loads
        inside lower-chord panels: array (k, nb, 4)."""
        states = []
        for index, chord in enumerate(self.chords):
            first, lengths = chord["first"], chord["lengths"]
            ne = len(lengths)
            ids = first + np.arange(ne)
            ui, uj = u[3 * ids], u[3 * ids + 3]
            axial = chord["ea"] * (uj - ui).T / lengths[None, :]
            ue = np.stack([u[3 * ids + 1], u[3 * ids + 2], u[3 * ids + 4], u[3 * ids + 5]], axis=2)   # (ne, k, 4)
            ends = np.einsum("eab,ekb->kea", chord["kb"], ue)
            if permanent:
                ends = ends + chord["fixed"][None, :, :]
            if extra is not None and index == 0:
                ends = ends + extra
            states.append({"N": axial, "w": ue[:, :, 0].T, "t": ue[:, :, 1].T,
                           "M": -ends[:, :, 1], "S": ends[:, :, 0]})
        return states

    # ---------------------------------------------------------- point forces
    def _locate(self, a):
        a = np.clip(np.asarray(a, dtype=float), 0.0, self.span)
        e = np.clip(np.searchsorted(self.xb, a, side="right") - 1, 0, self.nb - 1)
        t = np.clip(a - self.xb[e], 0.0, self.lb[e])
        return e, t

    def point_response(self, loads):
        """Superpose downward point forces using the exact clamped panel end actions.

        Each pair supplies an array of stations (one per case) and a force in N.
        """
        k = len(loads[0][0])
        rhs = np.zeros((self.nd, k))
        extra = np.zeros((k, self.nb, 4))
        cuts = []
        for stations, share in loads:
            e, t = self._locate(stations)
            l = self.lb[e]
            w_end = share * (l - t) ** 3 / 6.0
            t_end = share * (l - t) ** 2 / 2.0
            m0, s0 = _clamped_start(l, w_end, t_end)
            react = np.stack([s0, -m0, -(s0 - share), m0 + s0 * l - share * (l - t)], axis=1)
            cols = np.arange(k)
            for j, dof in enumerate((3 * e + 1, 3 * e + 2, 3 * e + 4, 3 * e + 5)):
                np.add.at(rhs, (dof, cols), -react[:, j])
            np.add.at(extra, (cols, e), react)
            cuts.append((e, t, share))
        u = self._solve(rhs)
        return u, self._chord_states(u, permanent=False, extra=extra), cuts

    # --------------------------------------------------------------- rating
    def envelope(self):
        if not hasattr(self, "_envelope"):
            from suspended_load import Envelope
            self._envelope = Envelope(self, SUSPENDED_FORCE / 1000., CLEARANCE, WHEELBASE, HEAVY_SHARE,
                                      STRESS_LIMIT, E_STEEL)
        return self._envelope

    def rating(self, stations):
        """Carrier capacity in N, safe for both orientations and every admissible suspended unit."""
        s = np.atleast_1d(np.asarray(stations, dtype=float))
        out = np.empty(len(s))
        envelope = self.envelope()
        # Bound transient array sizes without restricting the number of query stations.
        for part in np.array_split(np.arange(len(s)), max(1, math.ceil(len(s)/1024))):
            out[part] = 1000. * envelope.evaluate(s[part])[0]
        return out

    def state(self, station, heavy_at, unit_station):
        """Scalar demands for the explicitly located three point forces."""
        low, high = (HEAVY_SHARE, LIGHT_SHARE) if heavy_at == "LOW" else (LIGHT_SHARE, HEAVY_SHARE)
        s = np.array([station])
        u, live, wheels = self.point_response(((s, low), (s + WHEELBASE, high)))
        uh, fixed, unit = self.point_response(((np.array([unit_station]), SUSPENDED_FORCE),))
        perm = [{key: p[key] + h[key] for key in p} for p, h in zip(self.perm, fixed)]
        return _State(self, u[:, 0], live, wheels, self.web_perm + self.web_rows @ uh[:, 0], perm, unit[0])


class _State:
    def __init__(self, truss, u, live, wheels, web, perm, unit):
        self.truss, self.u, self.live, self.wheels = truss, u, live, wheels
        self.web, self.perm, self.unit = web, perm, unit

    def web_force(self, index, load):
        """Axial force, positive in tension, with steel, unit and carrier acting together."""
        return float(self.web[index] + load * (self.truss.web_rows[index] @ self.u))

    def _panel(self, chord_index, x, side):
        chord = self.truss.chords[chord_index]
        x0, lengths = chord["x0"], chord["lengths"]
        ends = np.concatenate([x0, [x0[-1] + lengths[-1]]])
        x = min(max(x, ends[0]), ends[-1])
        e = int(np.searchsorted(ends, x, side="right" if side > 0 else "left")) - 1
        e = min(max(e, 0), len(lengths) - 1)
        return chord, e, min(max(x - x0[e], 0.0), lengths[e])

    def _fields(self, chord_index, x, load, side):
        chord, e, t = self._panel(chord_index, x, side)
        ei, q = chord["ei"], chord["q"]
        pp, lp = self.perm[chord_index], self.live[chord_index]
        get = lambda st, key: float(st[key][0, e])
        mp = get(pp, "M") + get(pp, "S") * t - 0.5 * q * t * t
        wp = get(pp, "w") + get(pp, "t") * t + (get(pp, "M") * t * t / 2 + get(pp, "S") * t ** 3 / 6 - q * t ** 4 / 24) / ei
        ml = get(lp, "M") + get(lp, "S") * t
        wl = get(lp, "w") + get(lp, "t") * t + (get(lp, "M") * t * t / 2 + get(lp, "S") * t ** 3 / 6) / ei
        if chord_index == 0:
            ue, ut, weight = self.unit
            if int(ue[0]) == e:
                arm = max(t - float(ut[0]), 0.)
                mp -= weight * arm
                wp -= weight * arm ** 3 / (6. * ei)
            for we, wt, share in self.wheels:
                if int(we[0]) == e:
                    arm = max(t - float(wt[0]), 0.0)
                    ml -= share * arm
                    wl -= share * arm ** 3 / (6.0 * ei)
        axial = get(pp, "N") + load * get(lp, "N")
        return axial, mp + load * ml, wp + load * wl

    def chord_stress(self, chord_index, x, load, side):
        """Largest absolute fibre stress of a chord section at abscissa x."""
        chord = self.truss.chords[chord_index]
        member = chord["member"]
        axial, moment, _ = self._fields(chord_index, x, load, side)
        bend = moment * member["c"] / member["I"]
        return max(abs(axial / member["A"] + bend), abs(axial / member["A"] - bend))

    def sag(self, x, load):
        """Downward displacement of the lower chord axis at station x."""
        return -self._fields(0, x, load, 1)[2]


def load_trusses(path):
    """Read every TRUSS assembly of an IFC file."""
    import ifcopenshell
    model = ifcopenshell.open(path)
    to_mm = _length_factor(model)
    angle = _angle_factor(model)
    trusses = []
    cache = {}
    for assembly in model.by_type("IfcElementAssembly"):
        if assembly.PredefinedType != "TRUSS":
            continue
        members = []
        for rel in assembly.IsDecomposedBy:
            for obj in rel.RelatedObjects:
                if not obj.is_a("IfcMember"):
                    continue
                world = _object_placement(obj.ObjectPlacement)
                p1, p2 = _axis_ends(obj, world, to_mm)
                length = float(np.linalg.norm(p2 - p1))
                area, tensor, outline, centre = _body_section(obj, world, 0.5 * (p1 + p2), (p2 - p1) / length, to_mm, angle, cache)
                members.append({"guid": obj.GlobalId, "p1": p1, "p2": p2, "length": length, "A": float(area),
                                "tensor": tensor, "outline": outline, "centre": centre})
        trusses.append(Truss(assembly.GlobalId, members))
    return trusses
