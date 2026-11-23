#!/bin/bash
set -euo pipefail

mkdir -p /app/out
rm -f -- /app/out/rate_trusses.py /app/out/truss_ratings.json
cat > /app/out/rate_trusses.py <<'SPICA_RATING_PROGRAM'
#!/opt/conda/bin/python
"""Compute a continuously bounded two-wheel carrier rating, including a movable service unit.

Lengths are in mm. Structural forces are in N. Influence functions are for a
downward 1 kN load and positive utilization always means closer to a limit.
"""
import argparse
import io
import json
import math
from dataclasses import dataclass
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom as geom
from ifcopenshell import ifcopenshell_wrapper as ifc_kernel
import numpy as np
from OCP.Bnd import Bnd_Box
from OCP.BRep import BRep_Builder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.GProp import GProp_GProps
from OCP.TopoDS import TopoDS_Shape
from OCP.gp import gp_Dir, gp_Pln, gp_Pnt
from numba import njit
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

E = 200000.0
WEIGHT = 0.0000770
STRESS_LIMIT = 160.0
SPACING = 1100.0
JOINT_TOL = 0.05
UNIT_KN = 40.0
CLEARANCE = 400.0
KNOT_FACTOR = 0.9985
CHART_FLOOR = 0.992


@dataclass
class Member:
    guid: str
    ends: np.ndarray
    area: float
    inertia: float
    fibres: tuple | None  # bottom and top ordinates relative to the section centroid

    def record(self):
        return dict(guid=self.guid, ends_mm=self.ends.tolist(), A_mm2=self.area,
                    I_mm4=self.inertia, c_mm=None if self.fibres is None
                    else max(abs(z) for z in self.fibres))


class PlacedGeometry:
    def __init__(self):
        self.axis_settings = geom.settings()
        self.axis_settings.set("use-world-coords", True)
        self.axis_settings.set("dimensionality", ifc_kernel.CURVES)
        self.body_settings = geom.settings()
        self.body_settings.set("use-world-coords", True)
        self.body_settings.set("iterator-output", ifc_kernel.SERIALIZED)

    def member(self, product):
        reps = {r.RepresentationIdentifier.lower(): r
                for r in product.Representation.Representations}
        axis = geom.create_shape(self.axis_settings, product, reps["axis"])
        points = np.asarray(axis.geometry.verts).reshape(-1, 3) * 1000.0
        component = np.argmax(np.ptp(points, axis=0))
        ends = points[[points[:, component].argmin(), points[:, component].argmax()]]
        delta = ends[1] - ends[0]
        tangent = delta / np.linalg.norm(delta)

        # IfcOpenShell resolves mapped items, product placements and units.
        # Its serialized B-rep is in world metres. The common with an infinite
        # planar face is the actual steel section, including native conics.
        element = geom.create_shape(self.body_settings, product, reps["body"])
        body = TopoDS_Shape()
        BRepTools.Read_s(body, io.BytesIO(element.geometry.brep_data.encode()), BRep_Builder())
        point = gp_Pnt(*(ends.mean(axis=0) / 1000.0))
        plane = gp_Pln(point, gp_Dir(*tangent))
        section = BRepAlgoAPI_Common(body, BRepBuilderAPI_MakeFace(plane).Face())
        if not section.IsDone():
            raise ValueError(f"Body section failed for {product.GlobalId}")
        props = GProp_GProps(point)
        BRepGProp.SurfaceProperties_s(section.Shape(), props, 1e-10)
        area = props.Mass() * 1e6
        if area <= 0:
            raise ValueError(f"Empty Body at Axis midlength for {product.GlobalId}")
        matrix = props.MatrixOfInertia()  # about the section centroid, in world axes
        tensor = 1e12 * np.array([[matrix.Value(i, j) for j in (1, 2, 3)]
                                  for i in (1, 2, 3)])
        if abs(delta[2]) <= JOINT_TOL:
            normal = np.cross(tangent, [0.0, 0.0, 1.0])
            normal /= np.linalg.norm(normal)
            inertia = normal @ tensor @ normal
            box = Bnd_Box()
            BRepBndLib.AddOptimal_s(section.Shape(), box, False, False)
            zmin, zmax = box.Get()[2::3]
            zc = props.CentreOfMass().Z()
            fibres = ((zmin - zc) * 1000.0, (zmax - zc) * 1000.0)
        else:
            # For a plane section the normal principal moment is the polar sum.
            inertia = np.linalg.eigvalsh(tensor)[0]
            fibres = None
        return Member(product.GlobalId, ends, area, float(inertia), fibres)


class Frame:
    def __init__(self, members):
        self.chords = sorted((m for m in members if m.fibres is not None),
                             key=lambda m: m.ends[:, 2].mean())
        if len(self.chords) != 2:
            raise ValueError("Expected two horizontal chords")
        self.webs = [m for m in members if m.fibres is None]
        a, b = self.chords[0].ends
        roundoff = 8*np.finfo(float).eps*max(np.max(np.abs([a, b])), np.linalg.norm(b-a))
        component = 1 if abs(a[0] - b[0]) <= 0.1 + roundoff else 0
        start, end = (a, b) if a[component] < b[component] else (b, a)
        self.span = np.linalg.norm(end - start)
        along = (end - start) / self.span
        self.sag_limit = min(self.span / 600.0, 20.0)
        web_ends = np.array([w.ends for w in self.webs])
        web_delta = web_ends[:, 1] - web_ends[:, 0]
        posts = np.linalg.norm(web_delta[:, :2], axis=1) <= 2*JOINT_TOL + roundoff
        levels = np.array([m.ends[:, 2].mean() for m in self.chords])
        tiers = np.abs(web_ends[:, :, 2, None] - levels).argmin(axis=2)
        joints = np.empty(tiers.shape, dtype=int)
        positions, left_nodes, beam_tiers, beam_lengths = [], [], [], []
        for tier, chord in enumerate(self.chords):
            a, b = chord.ends
            if (b - a) @ along < 0:
                a, b = b, a
            length = np.linalg.norm(b - a)
            direction = (b - a) / length
            mask = tiers == tier
            # Every panel boundary has a post. Project its end onto this chord
            # and attach coincident diagonal ends to the same analytical node.
            post_ends = web_ends[mask & posts[:, None]]
            stations = np.sort((post_ends - a) @ direction)
            stations[[0, -1]] = [0.0, length]
            projections = (web_ends[mask] - a) @ direction
            first = len(positions)
            nearest = np.abs(projections[:, None] - stations).argmin(axis=1)
            joints[mask] = first + nearest
            positions.extend(a + stations[:, None] * direction)
            n = len(stations) - 1
            left_nodes.extend(range(first, first + n))
            beam_tiers.extend([tier] * n)
            beam_lengths.extend(np.diff(stations))
            if tier == 0:
                self.stations = stations
                self.n_lower = n

        self.positions = np.array(positions)
        self.left = np.array(left_nodes)
        self.tier = np.array(beam_tiers)
        self.length = np.array(beam_lengths)
        self.area = np.array([self.chords[t].area for t in self.tier])
        self.inertia = np.array([self.chords[t].inertia for t in self.tier])
        self.fibres = np.array([self.chords[t].fibres for t in self.tier])
        self.joints = joints

        ndof = 3 * len(positions)
        stiffness = np.zeros((ndof, ndof))
        dead = np.zeros(ndof)
        for i, (j, length) in enumerate(zip(self.left, self.length)):
            axial = np.array([3*j, 3*(j+1)])
            bend = np.array([3*j+1, 3*j+2, 3*j+4, 3*j+5])
            stiffness[np.ix_(axial, axial)] += E*self.area[i]/length * np.array([[1, -1], [-1, 1]])
            l2 = length**2
            kb = np.array([[12, 6*length, -12, 6*length],
                           [6*length, 4*l2, -6*length, 2*l2],
                           [-12, -6*length, 12, -6*length],
                           [6*length, 2*l2, -6*length, 4*l2]])
            stiffness[np.ix_(bend, bend)] += E*self.inertia[i]/length**3 * kb
            q = -WEIGHT*self.area[i]
            dead[bend] += q*np.array([length/2, l2/12, length/2, -l2/12])

        self.web_area = np.array([m.area for m in self.webs])
        web_length = np.linalg.norm(web_delta, axis=1)
        self.buckling = (math.pi**2 * E * np.array([w.inertia for w in self.webs])
                         / (2 * web_length**2))
        self.web_dofs = 3*joints[:, [0, 0, 1, 1]] + [0, 1, 0, 1]
        cx = (web_delta @ along) / web_length
        cz = web_delta[:, 2] / web_length
        extension = np.column_stack((-cx, -cz, cx, cz))
        self.web_recovery = (E*self.web_area/web_length)[:, None] * extension
        for i, dofs in enumerate(self.web_dofs):
            stiffness[np.ix_(dofs, dofs)] += np.outer(extension[i], self.web_recovery[i])
            dead[dofs[[1, 3]]] -= WEIGHT*self.web_area[i]*web_length[i]/2

        # One factorization serves gravity and every Hermite equivalent wheel load.
        lower_bend = (3*np.arange(self.n_lower + 1)[:, None] + [1, 2]).ravel()
        rhs = np.zeros((ndof, len(lower_bend) + 1))
        rhs[:, 0] = dead
        rhs[lower_bend, 1 + np.arange(len(lower_bend))] = 1.0
        free = np.ones(ndof, dtype=bool)
        free[[0, 1, 3*self.n_lower + 1]] = False
        solution = np.zeros_like(rhs)
        elastic = stiffness[np.ix_(free, free)]
        scale = 1/np.sqrt(np.diag(elastic))
        factor = cho_factor(scale[:, None]*elastic*scale)
        solution[free] = scale[:, None]*cho_solve(factor, scale[:, None]*rhs[free])
        self.dead_u = solution[:, 0]
        self.flexibility = solution[:, 1:]

    def crossing_stations(self):
        offsets = [0, SPACING, -CLEARANCE, CLEARANCE, SPACING-CLEARANCE, SPACING+CLEARANCE]
        last = self.span - SPACING
        events = sorted(x-d for x in self.stations for d in offsets if 0 < x-d < last)
        result = [0.0]
        for x in events:
            if x-result[-1] > 1e-9 and last-x > 1e-9:
                result.append(float(x))
        return result + [last]


def basis(n, t):
    t = np.asarray(t)
    return np.stack([math.comb(n, k)*t**k*(1-t)**(n-k) for k in range(n+1)], axis=-1)


def split_at(n, fraction):
    work = np.eye(n+1)
    left, right = [work[0]], [work[-1]]
    for i in range(n):
        work = (1-fraction)*work[:-1] + fraction*work[1:]
        left.append(work[0])
        right.append(work[-1])
    return np.array(left), np.array(right[::-1])


# The trivariate reserve has degrees <= (6, 6, 4) in carrier station,
# admissible unit position and section position in a load-free beam piece.
S6 = (1-np.cos(np.pi*np.arange(7)/6))/2
X4 = (1-np.cos(np.pi*np.arange(5)/4))/2
B6 = np.linalg.inv(basis(6, S6))
B4 = np.linalg.inv(basis(4, X4))
HALF6 = split_at(6, 0.5)
HALF4 = split_at(4, 0.5)
BINS = 8
EPS = 2e-13    # dimensionless utilization, roundoff in control coefficients


class Influences:
    """Section response for gravity and a downward 1 kN anywhere on the lower chord.

    dead[r] is in powers of the normalized section abscissa x. influence[r,p,k]
    is in powers of x and load abscissa y in lower panel p. The two polynomials
    k=0,1 give the extensions on either side of the point force in its own beam.
    """
    def __init__(self, f):
        self.frame = f
        n = f.n_lower
        columns = 1+4*n
        u = np.empty((len(f.dead_u), columns))
        u[:, 0] = f.dead_u
        for j, length in enumerate(f.length[:n]):
            hermite = np.array([[1, 0, -3, 2], [0, length, -2*length, length],
                                [0, 0, 3, -2], [0, 0, -length, length]])
            u[:, 1+4*j:5+4*j] = -1000*f.flexibility[:, 2*j:2*j+4] @ hermite
        self.checks, self.panel = [], []
        gravity, points = [], []

        def row(mode, kind, index, coeffs, jump=None):
            self.checks.append((mode, kind, index))
            self.panel.append(index if kind == "beam" and index < n else -1)
            gravity.append(coeffs[:, 0].copy())
            p = coeffs[:4, 1:].reshape(4, n, 4).transpose(1, 0, 2)
            both = np.repeat(p[:, None], 2, axis=1)
            if jump is not None and index < n:
                both[index, 1] += jump
            points.append(both)

        for i, (node, length) in enumerate(zip(f.left, f.length)):
            left = 3*node
            v0, rot0 = u[left+1], length*u[left+2]
            v1, rot1 = u[left+4], length*u[left+5]
            axial = (u[left+3]-u[left])/length
            disp = np.zeros((5, columns))
            disp[4, 0] = -WEIGHT*f.area[i]*length**4/(24*E*f.inertia[i])
            disp[0], disp[1] = v0, rot0
            disp[2] = 3*(v1-v0)-2*rot0-rot1+disp[4]
            disp[3] = -2*(v1-v0)+rot0+rot1-2*disp[4]
            # Integrate the clamped point-force particular solution. Behind
            # the load add k*(x-y)^3. This restores shear, moment, v and slope
            # throughout the beam, including between the interpolation nodes.
            k = -1000*length**3/(6*E*f.inertia[i])
            if i < n:
                cols = slice(1+4*i, 5+4*i)
                disp[2, cols] += 3*k*np.array([0, 1, -2, 1])
                disp[3, cols] -= k*np.array([1, 0, -3, 2])
            for fibre in f.fibres[i]:
                stress = np.zeros_like(disp)
                stress[:3] = (-E*fibre/(STRESS_LIMIT*length**2)
                              * np.array([2, 6, 12])[:, None] * disp[2:])
                stress[0] += E/STRESS_LIMIT*axial
                jump = np.zeros((4, 4))
                jump[1, 0] = -6*k*E*fibre/(STRESS_LIMIT*length**2)
                jump[0, 1] = -jump[1, 0]
                for sign in (1, -1):
                    row("STRESS", "beam", i, sign*stress, sign*jump)
            if i < n:
                jump = np.zeros((4, 4))
                for d in range(4):
                    jump[3-d, d] = -k/f.sag_limit*math.comb(3, d)*(-1)**d
                row("SAG", "beam", i, -disp/f.sag_limit, jump)

        for i, dofs in enumerate(f.web_dofs):
            force = f.web_recovery[i] @ u[dofs]
            allowable = STRESS_LIMIT*f.web_area[i]
            coeffs = np.zeros((5, columns))
            coeffs[0] = force/allowable
            row("STRESS", "web", i, coeffs)
            if f.buckling[i] < allowable:
                coeffs[0] = -force/f.buckling[i]
                row("WEB_BUCKLING", "web", i, coeffs)
            else:
                coeffs[0] = -force/allowable
                row("STRESS", "web", i, coeffs)
        self.panel = np.array(self.panel, dtype=np.int64)
        self.dead = np.ascontiguousarray(gravity)
        self.effect = np.ascontiguousarray(points)
        self.dead_max = (self.dead @ (np.vander(X4, 5, increasing=True).T @ B4.T)).max(axis=1)
        anchors = np.array([0., 0.5, 1.])
        self.dead_anchor = self.dead @ np.vander(anchors, 5, increasing=True).T
        self.upper, self.lower_anchor, self.unit_anchor = influence_bounds(self.effect, self.panel)
        self.anchor_s = np.r_[np.concatenate([f.stations[i]+f.length[i]*np.arange(BINS)/BINS
                                              for i in range(n)]), f.span]


@njit
def power_values(coeff, y):
    return coeff[:, 0] + y*(coeff[:, 1] + y*(coeff[:, 2] + y*coeff[:, 3]))


@njit(inline="always")
def cubic(coeff, x):
    return coeff[0]+x*(coeff[1]+x*(coeff[2]+x*coeff[3]))


@njit(inline="always")
def quartic(coeff, x):
    return coeff[0]+x*(coeff[1]+x*(coeff[2]+x*(coeff[3]+x*coeff[4])))


@njit
def influence_bounds(effect, panel):
    """Certified coarse bounds used solely to discard individually harmless checks."""
    nr, npanels = effect.shape[:2]
    upper = np.empty((nr, npanels, BINS))
    anchors = np.empty((nr, npanels*BINS+1, 3))
    lower = np.empty((nr, npanels, BINS, 3))
    for r in range(nr):
        for p in range(npanels):
            local = panel[r] == p
            for b in range(BINS):
                highest = -np.inf
                at_fixed_x = np.empty((7, 3))
                for side in range(1+int(local)):
                    values = np.empty((7, 5))
                    for j in range(7):
                        y = (b+S6[j])/BINS
                        coeff = power_values(effect[r, p, side], y)
                        for k in range(5):
                            x = X4[k]
                            if local:
                                x = y+(1-y)*x if side else y*x
                            values[j, k] = cubic(coeff, x)
                        for k in range(3):
                            x = k*0.5
                            # Each fixed anchor is a bin boundary, so this
                            # bin's entire polynomial uses one branch.
                            if side == int(local and x >= (b+0.5)/BINS):
                                at_fixed_x[j, k] = cubic(coeff, x)
                    highest = max(highest, np.max(B6 @ values @ B4.T))
                upper[r, p, b] = highest
                controls = B6 @ at_fixed_x
                for k in range(3):
                    lower[r, p, b, k] = np.min(controls[:, k])
                    anchors[r, p*BINS+b, k] = at_fixed_x[0, k]
                    if p == npanels-1 and b == BINS-1:
                        anchors[r, npanels*BINS, k] = at_fixed_x[-1, k]
    return upper, lower, anchors


def bin_range(lo, hi):
    first = max(0, min(BINS-1, math.floor(lo*BINS)))
    last = max(first, min(BINS-1, math.ceil(hi*BINS)-1))
    return slice(first, last+1)


@njit
def samples_on_pieces(tags, edges, dead, effect):
    """Sample exact polynomials on a product of three unit intervals."""
    out = np.empty((len(tags), 3, 7, 7, 5))
    for g in range(len(tags)):
        r, p0, p1, ph, xl, xr, k0, k1, kh = tags[g]
        e = edges[g]
        where = np.array([0., 1., 0., 0., 0.])
        for i in range(7):
            t = S6[i]
            where[2] = (1-t)*e[0]+t*e[1]
            where[3] = (1-t)*e[2]+t*e[3]
            lower = (1-t)*e[4]+t*e[5]
            upper = (1-t)*e[6]+t*e[7]
            c0 = power_values(effect[r, p0, k0], where[2])
            c1 = power_values(effect[r, p1, k1], where[3])
            for j in range(7):
                where[4] = (1-S6[j])*lower+S6[j]*upper
                ch = power_values(effect[r, ph, kh], where[4])
                for k in range(5):
                    x = (1-X4[k])*where[xl]+X4[k]*where[xr]
                    m0, m1 = cubic(c0, x), cubic(c1, x)
                    out[g, 0, i, j, k] = 1-quartic(dead[r], x)-UNIT_KN*cubic(ch, x)
                    out[g, 1, i, j, k] = 0.65*m0+0.35*m1
                    out[g, 2, i, j, k] = 0.35*m0+0.65*m1
    return out


class StationCell:
    """The ordering of the three loads and the panel/clearance edges is fixed."""
    def __init__(self, inf, a, b):
        self.inf, self.a, self.b = inf, a, b
        f = inf.frame
        ends = np.array([a, b])
        p0, p1 = np.searchsorted(f.stations, (a+b)/2+np.array([0., SPACING]), side="right")-1
        wheel0 = (ends-f.stations[p0])/f.length[p0]
        wheel1 = (ends+SPACING-f.stations[p1])/f.length[p1]
        bins0, bins1 = bin_range(*wheel0), bin_range(*wheel1)
        w0_lo = inf.lower_anchor[:, p0, bins0].min(axis=1)
        w1_lo = inf.lower_anchor[:, p1, bins1].min(axis=1)
        m_lo = np.maximum(.65*w0_lo+.35*w1_lo, .35*w0_lo+.65*w1_lo)
        # A section and a unit station common to the whole cell give an
        # attained upper capacity. This ceiling is only a pruning device.
        hs = inf.anchor_s
        common = (hs <= a-CLEARANCE) | (hs >= b+SPACING+CLEARANCE)
        common |= (hs >= b+CLEARANCE) & (hs <= a+SPACING-CLEARANCE)
        fixed_demand = inf.dead_anchor + UNIT_KN*np.max(inf.unit_anchor[:, common, :], axis=1)
        upper_load = np.full_like(m_lo, np.inf)
        np.divide(1-fixed_demand+EPS, m_lo-EPS, where=m_lo>EPS, out=upper_load)
        self.ceiling = max(0., float(upper_load.min()))*(1+1e-9)+1e-10

        w0_hi = inf.upper[:, p0, bins0].max(axis=1)
        w1_hi = inf.upper[:, p1, bins1].max(axis=1)
        m_hi = np.maximum(.65*w0_hi+.35*w1_hi, .35*w0_hi+.65*w1_hi)
        finite_ceiling = np.isfinite(self.ceiling)
        always = (inf.dead_max + self.ceiling*np.maximum(m_hi, 0.) if finite_ceiling
                  else np.full_like(m_hi, np.inf))
        pieces, edges = [], []
        windows = [(np.zeros(2), ends-CLEARANCE, (4, 2, 3)),
                   (ends+CLEARANCE, ends+SPACING-CLEARANCE, (2, 4, 3)),
                   (ends+SPACING+CLEARANCE, np.full(2, f.span), (2, 3, 4))]
        for lo, hi, order in windows:
            for ph in range(f.n_lower):
                start, length = f.stations[ph], f.length[ph]
                bottom = np.maximum(0., (lo-start)/length)
                top = np.minimum(1., (hi-start)/length)
                if top.mean() <= bottom.mean():
                    continue
                bins = bin_range(min(bottom), max(top))
                possible = always + UNIT_KN*inf.upper[:, ph, bins].max(axis=1)
                for r in np.flatnonzero(possible >= 1-3e-11):
                    locations = (p0, p1, ph)
                    cuts = [0]+[i for i in order if locations[i-2] == inf.panel[r]]+[1]
                    sides = [0, 0, 0]
                    for left, right in zip(cuts, cuts[1:]):
                        if left > 1:
                            sides[left-2] = 1
                        pieces.append((r, p0, p1, ph, left, right, *sides))
                        edges.append((*wheel0, *wheel1, *bottom, *top))

        tags = np.array(pieces, dtype=np.int64)
        positions = np.array(edges)
        nets, selected, orientations = [], [], []
        for start in range(0, len(tags), 96):
            stop = min(len(tags), start+96)
            values = samples_on_pieces(tags[start:stop], positions[start:stop], inf.dead, inf.effect)
            coeff = values @ B4.T
            coeff = B6 @ coeff
            coeff = (B6 @ coeff.reshape(-1, 7, 35)).reshape(coeff.shape)
            h = coeff[:, 0]
            for orientation in (0, 1):
                m = coeff[:, 1+orientation]
                keep = ((h-self.ceiling*np.maximum(m, 0.)).min(axis=(1, 2, 3)) <= 3e-11
                        if finite_ceiling else np.ones(len(h), dtype=bool))
                if np.any(keep):
                    nets.append(np.stack((h[keep], m[keep]), axis=1))
                    selected.extend(start+np.flatnonzero(keep))
                    orientations.extend([orientation]*int(keep.sum()))
        if not nets:
            raise ArithmeticError("No potentially active section")
        self.tags = tags[selected]
        self.edges = positions[selected]
        self.orient = np.array(orientations)
        # Axes: station, oriented section piece, reserve/influence, unit, section.
        self.net = np.ascontiguousarray(np.concatenate(nets).transpose(2, 0, 1, 3, 4))

    def witness(self, s, lim):
        r, p0, p1, ph, xl, xr, *_ = self.tags[lim.group]
        e = self.edges[lim.group]
        u = (s-self.a)/(self.b-self.a)
        wheel0, wheel1, bot, top = (1-u)*e[::2]+u*e[1::2]
        h = (1-lim.v)*bot+lim.v*top
        positions = [0., 1., wheel0, wheel1, h]
        x = (1-lim.z)*positions[xl]+lim.z*positions[xr]
        f = self.inf.frame
        mode, kind, i = self.inf.checks[r]
        if kind == "beam":
            member = f.chords[f.tier[i]]
            j = f.left[i]
            xyz = (1-x)*f.positions[j]+x*f.positions[j+1]
        else:
            member = f.webs[i]
            xyz = member.ends.mean(axis=0)
        return dict(s_mm=float(s), unit_s_mm=float(f.stations[ph]+h*f.length[ph]),
                    heavy_at=("LOW", "HIGH")[self.orient[lim.group]], mode=mode,
                    member=member.guid, xyz_mm=xyz.tolist())


@dataclass
class Limit:
    load: float
    group: int
    v: float             # fraction across the admissible part of the unit panel
    z: float             # fraction along the load-free section piece


@njit
def improve_ratio(h, m, g, vl, vr, zl, zr, best, rel):
    """de Casteljau branch and bound for the attained minimum of H/M, M>0."""
    for i in (0, 6):
        for j in (0, 4):
            if m[i, j] > 0. and h[i, j] < best[0]*m[i, j]:
                best[0] = h[i, j]/m[i, j]
                best[1] = g
                best[2] = vl if i == 0 else vr
                best[3] = zl if j == 0 else zr
    diff = h-(1-rel)*best[0]*m
    if np.min(diff) >= -EPS:
        return
    dh = 6*np.max(np.abs(diff[1:]-diff[:-1]))
    dx = 4*np.max(np.abs(diff[:, 1:]-diff[:, :-1]))
    if dh > dx:
        mid = (vl+vr)/2
        improve_ratio(HALF6[0]@h, HALF6[0]@m, g, vl, mid, zl, zr, best, rel)
        improve_ratio(HALF6[1]@h, HALF6[1]@m, g, mid, vr, zl, zr, best, rel)
    else:
        mid = (zl+zr)/2
        improve_ratio(h@HALF4[0].T, m@HALF4[0].T, g, vl, vr, zl, mid, best, rel)
        improve_ratio(h@HALF4[1].T, m@HALF4[1].T, g, vl, vr, mid, zr, best, rel)


@njit
def ratio_minimum(net, rel):
    best = np.array([np.inf, 0., 0., 0.])
    for g in range(net.shape[0]):
        for i in (0, 6):
            for j in (0, 4):
                h, m = net[g, 0, i, j], net[g, 1, i, j]
                if m > 0 and h/m < best[0]:
                    best[0], best[1], best[2], best[3] = h/m, g, i/6., j/4.
    if not np.isfinite(best[0]):
        # An interior-only live effect also supplies an attained starting bound.
        for g in range(net.shape[0]):
            v, z = negative_point(-net[g, 1], 0., 1., 0., 1.)
            if v >= 0:
                bv = CHOOSE6*v**np.arange(7)*(1-v)**np.arange(6, -1, -1)
                bz = CHOOSE4*z**np.arange(5)*(1-z)**np.arange(4, -1, -1)
                h = bv @ net[g, 0] @ bz
                m = bv @ net[g, 1] @ bz
                if m > 0 and h/m < best[0]:
                    best[0], best[1], best[2], best[3] = h/m, g, v, z
    if np.isfinite(best[0]):
        for g in range(net.shape[0]):
            improve_ratio(net[g, 0], net[g, 1], g, 0., 1., 0., 1., best, rel)
    return best


CHOOSE6 = np.array([math.comb(6, i) for i in range(7)])
CHOOSE4 = np.array([math.comb(4, i) for i in range(5)])


def station_limit(section):
    load, group, v, z = ratio_minimum(section, 2e-7)
    if not np.isfinite(load):
        raise ArithmeticError("No attained stationary rating")
    return Limit(float(load), int(group), float(v), float(z))


@njit
def negative_point(row, vl, vr, zl, zr):
    """Return a point of negative reserve, or certify the whole rectangle."""
    if np.min(row) >= -EPS:
        return -1., -1.
    for i in (0, 6):
        for j in (0, 4):
            if row[i, j] < -EPS:
                return (vl if i == 0 else vr), (zl if j == 0 else zr)
    dv = 6*np.max(np.abs(row[1:]-row[:-1]))
    dz = 4*np.max(np.abs(row[:, 1:]-row[:, :-1]))
    if dv > dz:
        mid = (vl+vr)/2
        v, z = negative_point(HALF6[0] @ row, vl, mid, zl, zr)
        if v < 0:
            v, z = negative_point(HALF6[1] @ row, mid, vr, zl, zr)
    else:
        mid = (zl+zr)/2
        v, z = negative_point(row @ HALF4[0].T, vl, vr, zl, mid)
        if v < 0:
            v, z = negative_point(row @ HALF4[1].T, vl, vr, mid, zr)
    return v, z


@njit
def first_negative(rows):
    for g in range(len(rows)):
        v, z = negative_point(rows[g], 0., 1., 0., 1.)
        if v >= 0:
            return g, v, z
    return -1, -1., -1.


def nonnegative(reserve):
    return first_negative(np.ascontiguousarray(reserve.reshape(-1, 7, 5)))[0] < 0


def split_station(net, fraction=0.5):
    matrices = HALF6 if fraction == 0.5 else split_at(6, fraction)
    flat = net.reshape(7, -1)
    return tuple((m @ flat).reshape(net.shape) for m in matrices)


def line_reserve(h, m, low, high):
    """Degree-seven station coefficients of H minus a linear load times M."""
    q = h-low*m
    r = h-high*m
    w = (np.arange(1, 7)/7).reshape((6,)+(1,)*(h.ndim-1))
    return np.concatenate((q[:1], (1-w)*q[1:]+w*r[:-1], r[-1:]))


def polish(cell, s, lim, fixed_station=False):
    """Improve an attained minimum inside its own closed polynomial piece."""
    h, m = np.moveaxis(cell.net[:, lim.group], 1, 0)
    u0 = (s-cell.a)/(cell.b-cell.a)
    def evaluate(coords):
        u, v, z = coords
        bs = [basis(6, u), basis(6, v), basis(4, z)]
        grad = []
        def dot(field, a, b, c):
            return a @ (b @ field @ c)
        hm = np.array([dot(d, *bs) for d in (h, m)])
        if hm[1] <= 0:
            return 1e20, np.zeros(3)
        value = hm[0]/hm[1]
        residual = (h-value*m)/hm[1]
        for j, (t, degree) in enumerate(zip(coords, (6, 6, 4))):
            b = bs.copy()
            b[j] = basis(degree-1, t)
            grad.append(degree*dot(np.diff(residual, axis=j), *b))
        return value, np.array(grad)
    bounds = [(u0, u0) if fixed_station else (0., 1.), (0., 1.), (0., 1.)]
    fit = minimize(evaluate, [u0, lim.v, lim.z], jac=True, method="L-BFGS-B", bounds=bounds,
                   options=dict(ftol=1e-15, gtol=1e-10, maxiter=100))
    if np.isfinite(fit.fun) and fit.fun < lim.load:
        u, v, z = fit.x
        location = s if fixed_station else (1-u)*cell.a+u*cell.b
        return location, Limit(float(fit.fun), lim.group, v, z)
    return s, lim


def minimum_in_cell(cell, previous, first, last):
    """Locate a continuous minimizer capable of improving the best preceding cell."""
    target = previous
    current = None
    def consider(s, lim):
        nonlocal target, current
        if lim.load < target:
            s, lim = polish(cell, s, lim)
            target, current = lim.load, (s, lim)
    consider(cell.a, first)
    consider(cell.b, last)

    def search(a, b, net):
        if nonnegative(net[:, :, 0]-(1-2e-8)*target*net[:, :, 1]):
            return
        left, right = split_station(net)
        mid = (a+b)/2
        consider(mid, station_limit(left[-1]))
        if b-a > 0.12:
            search(a, mid, left)
            search(mid, b, right)
    search(cell.a, cell.b, cell.net)
    return current


def draw_chart(cell, a, b, net, left_limit, right_limit, knots):
    """Prove both sides of the load band over every closed chart interval."""
    low, high = KNOT_FACTOR*left_limit.load, KNOT_FACTOR*right_limit.load
    left, right = split_station(net)
    mid = (a+b)/2
    middle = station_limit(left[-1])
    # Along a fixed (unit fraction, section fraction) path, the midstation
    # configuration remains admissible everywhere. If C/CHART_FLOOR activates
    # that check throughout the interval, the entire chart line is high enough.
    h, m = (basis(6, middle.v) @ net[:, middle.group] @ basis(4, middle.z)).T
    full = np.min(m) > 0 and np.max(line_reserve(h, m, low/CHART_FLOOR, high/CHART_FLOOR)) <= EPS
    if full and nonnegative(line_reserve(net[:, :, 0], net[:, :, 1], low, high)):
        for s, load, lim in ((a, low, left_limit), (b, high, right_limit)):
            if s not in knots or load < knots[s][0]:
                knots[s] = (load, cell.witness(s, lim))
        return
    if mid == a or mid == b:
        raise ArithmeticError("Cannot resolve carrier station interval")
    draw_chart(cell, a, mid, left, left_limit, middle, knots)
    draw_chart(cell, mid, b, right, middle, right_limit, knots)


def rate(members):
    frame = Frame(members)
    inf = Influences(frame)
    breaks = frame.crossing_stations()
    knots = {}
    global_best = np.inf
    for a, b in zip(breaks, breaks[1:]):
        cell = StationCell(inf, a, b)
        if a == 0:
            group, v, z = first_negative(np.ascontiguousarray(cell.net[0, :, 0]))
            if group >= 0:
                # Under the specified reissue contract one overloaded station
                # implies an admissible permanent overload at every station.
                return [[0., 0.], [breaks[-1], 0.]], cell.witness(0., Limit(0., group, v, z))
        first = station_limit(cell.net[0])
        last = station_limit(cell.net[-1])
        result = minimum_in_cell(cell, global_best, first, last)
        if result is not None:
            s, lim = result
            global_best = lim.load
            if a+1e-8 < s < b-1e-8:
                left, right = split_station(cell.net, (s-a)/(b-a))
                draw_chart(cell, a, s, left, first, lim, knots)
                draw_chart(cell, s, b, right, lim, last, knots)
                continue
            if s-a <= 1e-8:
                _, first = polish(cell, a, first, fixed_station=True)
            else:
                _, last = polish(cell, b, last, fixed_station=True)
        draw_chart(cell, a, b, cell.net, first, last, knots)
    s = min(knots, key=lambda t: knots[t][0])
    return [[t, knots[t][0]] for t in sorted(knots)], knots[s][1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    model = ifcopenshell.open(str(args.structure))
    cad = PlacedGeometry()
    result = []
    for assembly in model.by_type("IfcElementAssembly"):
        if assembly.PredefinedType == "TRUSS":
            members = [cad.member(obj) for rel in assembly.IsDecomposedBy
                       if rel.is_a("IfcRelAggregates") for obj in rel.RelatedObjects
                       if obj.is_a("IfcMember")]
            curve, witness = rate(members)
            result.append(dict(guid=assembly.GlobalId, members=[m.record() for m in members],
                               curve=curve, witness=witness))
    target = args.out / "truss_ratings.json"
    target.write_text(json.dumps({"trusses": result}, allow_nan=False, indent=2)+"\n",
                      encoding="utf-8")


if __name__ == "__main__":
    with threadpool_limits(limits=1, user_api="blas"):
        main()
SPICA_RATING_PROGRAM
/opt/conda/bin/python /app/out/rate_trusses.py --structure /app/inputs/model_7_structure.ifc --out /app/out
