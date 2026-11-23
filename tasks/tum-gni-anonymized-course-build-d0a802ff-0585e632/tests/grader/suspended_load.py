"""Envelope of a carrier and a separately positioned load on the lower chord.

For a fixed observation section, the influence line for a downward point load
is cubic in each load panel. The section itself is another break in that line.
Its worst admissible ordinate is at a stationary point, a break or a clearance
boundary. The wheel orientation maximising the same signed response is then
known, so the remaining capacity is a one-variable linear inequality.

The spatial search retains an attained upper bound on the permissible load and
certifies a lower bound over all sections. For stress, start with the cuts at
joints and at all three loads. Between cuts the signed stress is quadratic due
only to self-weight. At the cut under the suspended unit both its station and
the section move together. The permanent stress on this diagonal is quartic.
For sag, any load satisfying the stress checks has |w''| <= stress_limit/(E*c).
These curvature bounds enclose the departure from endpoint interpolation by
max |f''| times interval length squared over eight. Only undecided section
intervals are bisected. No grid of suspended-unit positions is used.
"""
import math
from typing import NamedTuple

import numpy as np
from numba import njit

# Numerical enclosure of the reference capacity, inside the model comparison allowance.
RELATIVE_GAP = 2.e-6


class Influence(NamedTuple):
    coeff: np.ndarray
    panel: int
    left: np.ndarray
    section: float
    peaks: np.ndarray
    heights: np.ndarray


@njit
def value(coeff, u):
    return ((coeff[3] * u + coeff[2]) * u + coeff[1]) * u + coeff[0]


@njit
def at(row, panel, u, station):
    if panel == row.panel and station < row.section:
        return value(row.left, u)
    return value(row.coeff[panel], u)


@njit
def lookup(row, joints, z):
    j = min(max(np.searchsorted(joints, z, side="right") - 1, 0), len(joints) - 2)
    return at(row, j, (z - joints[j]) / (joints[j + 1] - joints[j]), z)


@njit
def prepare_line(coeff, panel, left, section, joints):
    # Include both endpoints and all real derivative roots in each closed piece.
    pos = np.empty(4 * (len(coeff) + 1))
    mag = np.empty_like(pos)
    n = 0
    for j in range(len(coeff)):
        length = joints[j + 1] - joints[j]
        split = (section - joints[j]) / length if j == panel else 0.
        for side in range(2 if j == panel else 1):
            c = left if j == panel and side == 0 else coeff[j]
            a = split if side == 1 else 0.
            b = split if j == panel and side == 0 else 1.
            qa, qb, qc = 3. * c[3], 2. * c[2], c[1]
            roots = [a, b]
            if qa == 0.:
                if qb != 0.:
                    roots.append(-qc / qb)
            else:
                disc = qb * qb - 4. * qa * qc
                if disc >= 0.:
                    t = -.5 * (qb + math.copysign(math.sqrt(disc), qb))
                    if t != 0.:
                        roots.append(t / qa)
                        roots.append(qc / t)
            for u in roots:
                if a <= u <= b:
                    pos[n] = joints[j] + length * u
                    mag[n] = value(c, u)
                    n += 1
    order = np.argsort(-mag[:n])
    return Influence(coeff, panel, left, section, pos[:n][order], mag[:n][order])


@njit
def row_samples(row, selected, stations, panels, fractions, clearance, span):
    """Wheel influence ordinates and the attained worst feasible unit ordinate and station."""
    result = np.empty((len(selected), 4))
    for i, k in enumerate(selected):
        s, far = stations[k, 0], stations[k, 1]
        for j in range(2):
            result[i, j] = at(row, panels[k, j], fractions[k, j], stations[k, j])
        high, where = -np.inf, 0.
        for p in range(len(row.peaks)):
            z = row.peaks[p]
            if abs(z - s) >= clearance and abs(z - far) >= clearance:
                high, where = row.heights[p], z
                break
        if high < row.heights[0]:
            for j in range(2, 6):
                z = stations[k, j]
                if 0. <= z <= span:
                    y = at(row, panels[k, j], fractions[k, j], z)
                    if y > high:
                        high, where = y, z
        result[i, 2], result[i, 3] = high, where
    return result


@njit
def capacity(base, slope):
    if base > 1.:
        return 0.
    if slope > 0.:
        return (1. - base) / slope
    return np.inf


@njit
def wheel_rate(first, second, heavy):
    return heavy * max(first, second) + (1. - heavy) * min(first, second)


@njit
def observe(best, detail, k, code, x, h, first, second, base, heavy):
    bound = capacity(base, wheel_rate(first, second, heavy))
    if bound < best[k]:
        best[k] = bound
        detail[k, 0] = code
        detail[k, 1] = x
        detail[k, 2] = h
        detail[k, 3] = 0 if first >= second else 1


@njit
def endpoints(best, detail, ids, code, x, dead, samples, weight, heavy):
    for i, k in enumerate(ids):
        f, g, peak, h = samples[i]
        observe(best, detail, k, code, x, h, f, g, dead + weight * peak, heavy)


@njit
def tent(a, b, x, z):
    # Difference between the moment and its endpoint interpolant for a downward point force.
    if z <= a or z >= b:
        return 0.
    return (min(x, z) - a) * (b - max(x, z)) / (b - a)


@njit
def inspect_cell(ids, best, detail, stations, joints, code, bounds, dead, geometry, poly, rows, samples,
                 weight, clearance, heavy):
    """Certify the section interval at these stations, returning just those needing subdivision."""
    a, b = bounds
    da, db = dead
    beta, curvature, origin, panel_length, diagonal_curvature = geometry
    # beta = c/I per kN, divided by the stress limit, for positive moment on the lower chord.
    # curvature is a bound on the spatial second derivative between all physical cuts.
    length = b - a
    error = curvature * length * length / 8.
    first_row, last_row = rows
    sa, sb = samples
    keep = np.empty(len(ids), dtype=np.int64)
    count = 0
    for i, k in enumerate(ids):
        if best[k] == 0.:
            continue
        low = np.inf
        for end in range(2):
            rr = sa[i] if end == 0 else sb[i]
            base = (da if end == 0 else db) + weight * rr[2]
            low = min(low, capacity(base + error, wheel_rate(rr[0], rr[1], heavy)))
        if beta > 0.:
            s, far = stations[k, 0], stations[k, 1]
            # A cut at either wheel. The wheel field has its exact tent correction.
            # For the unit field endpoint maxima are an upper interpolant outside
            # the cell. Inside, its positive tent has height at most length / 4.
            # A cell narrower than the clearance cannot contain this wheel AND the unit.
            for j in range(2):
                x = stations[k, j]
                if not a < x < b:
                    continue
                v = (x - a) / length
                f = (1. - v) * sa[i, 0] + v * sb[i, 0] + beta * tent(a, b, x, s)
                g = (1. - v) * sa[i, 1] + v * sb[i, 1] + beta * tent(a, b, x, far)
                dx = (1. - v) * da + v * db + curvature * (x - a) * (b - x) / 2.
                above = (1. - v) * sa[i, 2] + v * sb[i, 2]
                if length > clearance:
                    above += beta * length / 4.
                low = min(low, capacity(dx + weight * above + error, wheel_rate(f, g, heavy)))
                # Both endpoint maximisers are feasible unit placements for this wheel section.
                ha, hb = sa[i, 3], sb[i, 3]
                if ha == hb:
                    actual = (1. - v) * sa[i, 2] + v * sb[i, 2] + beta * tent(a, b, x, ha)
                    observe(best, detail, k, code, x, ha, f, g, dx + weight * actual, heavy)
                else:
                    for h in (ha, hb):
                        actual = ((1. - v) * lookup(first_row, joints, h) + v * lookup(last_row, joints, h)
                                  + beta * tent(a, b, x, h))
                        observe(best, detail, k, code, x, h, f, g, dx + weight * actual, heavy)
            # Cut directly under the unit. On each feasible interval the wheel
            # fields are affine in section station and the other field is quartic.
            for left, right in ((a, min(b, s - clearance)), (max(a, s + clearance), min(b, far - clearance)),
                                (max(a, far + clearance), b)):
                if left > right:
                    continue
                de = diagonal_curvature * (right - left) ** 2 / 8.
                for x in (left, right):
                    v = (x - a) / length
                    f = (1. - v) * sa[i, 0] + v * sb[i, 0] + beta * tent(a, b, x, s)
                    g = (1. - v) * sa[i, 1] + v * sb[i, 1] + beta * tent(a, b, x, far)
                    u = (x - origin) / panel_length
                    base = (((poly[4] * u + poly[3]) * u + poly[2]) * u + poly[1]) * u + poly[0]
                    observe(best, detail, k, code, x, x, f, g, base, heavy)
                    low = min(low, capacity(base + de + error, wheel_rate(f, g, heavy)))
        if low < best[k] * (1. - RELATIVE_GAP):
            keep[count] = i
            count += 1
    return keep[:count]


class SectionField:
    """Signed, limit-normalised response in one chord panel."""

    def __init__(self, truss, coefficients, index, panel, axial, bending, sag, weight, stress, modulus):
        chord = truss.chords[index]
        self.joints = truss.xb
        self.origin = float(chord["x0"][panel])
        self.length = float(chord["lengths"][panel])
        self.local = panel if index == 0 else -1
        self.axial, self.bending, self.sag = axial, bending, sag
        self.ei, self.q = chord["ei"], chord["q"]
        self.perm = {key: val[0, panel] for key, val in truss.perm[index].items()}
        self.coef = {key: val[panel] for key, val in coefficients.items()}
        self.cache = {}
        self.diagonal = np.zeros(5)
        self.beta = 1000. * bending if self.local >= 0 and bending > 0. else 0.
        self.curvature = (stress / (modulus * chord["member"]["c"] * truss.sag_limit)
                          if sag else max(0., bending * self.q))
        self.diagonal_curvature = 0.
        if self.beta:
            # Evaluate the point-load start fields at z=x. The local lever arm then vanishes.
            a, L, p = self.coef, self.length, self.perm
            self.diagonal[:4] = weight * (axial * a["N"][panel] + bending * a["M"][panel])
            self.diagonal[1:] += weight * bending * L * a["S"][panel]
            self.diagonal[:3] += [axial * p["N"] + bending * p["M"], bending * L * p["S"],
                                  -bending * self.q * L**2 / 2.]
            c = self.diagonal
            self.diagonal_curvature = (2. * abs(c[2]) + 6. * abs(c[3]) + 12. * abs(c[4])) / L**2

    def row(self, x):
        if x in self.cache:
            return self.cache[x]
        t = x - self.origin
        a, p = self.coef, self.perm
        coeff = self.axial * a["N"] + self.bending * (a["M"] + t * a["S"])
        dead = self.axial * p["N"] + self.bending * (p["M"] + t * p["S"] - self.q * t*t / 2.)
        if self.sag:
            coeff = self.sag * (a["w"] + t * a["t"] + (t*t/2. * a["M"] + t**3/6. * a["S"]) / self.ei)
            dead = self.sag * (p["w"] + t * p["t"] + (t*t/2. * p["M"] + t**3/6. * p["S"]
                                                       - self.q * t**4/24.) / self.ei)
        left = np.zeros(4)
        if self.local >= 0:
            L = self.length
            left[:] = coeff[self.local]
            left[:2] -= 1000. * self.bending * np.array([t, -L])
            if self.sag:
                left -= 1000. * self.sag / (6. * self.ei) * np.array([t**3, -3.*t*t*L, 3.*t*L*L, -L**3])
        result = (prepare_line(np.ascontiguousarray(coeff), self.local, left, x, self.joints), dead)
        self.cache[x] = result
        return result


class Envelope:
    def __init__(self, truss, weight, clearance, wheelbase, heavy, stress, modulus):
        self.truss = truss
        self.weight, self.clearance, self.wheelbase, self.heavy = weight, clearance, wheelbase, heavy
        # Four interior evaluations recover every cubic exactly. Avoid assigning a point
        # at a joint to either adjacent panel while constructing the power coefficients.
        v = np.array([.125, .375, .625, .875])
        h = (truss.xb[:-1, None] + truss.lb[:, None] * v).ravel()
        u, states, _ = truss.point_response(((h, 1000.),))
        transform = np.linalg.inv(np.vander(v, 4, increasing=True)).T

        def power(values):
            return np.ascontiguousarray(values.reshape(len(values), truss.nb, 4) @ transform)

        self.web_lines = []
        self.fields = []
        self.labels = []
        webs = power(truss.web_rows @ u)
        for j, web in enumerate(truss.webs):
            strength = stress * web["A"]
            for sign in (1., -1.):
                limit = strength if sign > 0. else min(strength, truss.web_crit[j])
                mode = "WEB_BUCKLING" if limit < strength else "STRESS"
                self.labels.append((web, mode))
                row = prepare_line(webs[j] * (sign / limit), -1, np.zeros(4), 0., truss.xb)
                self.web_lines.append((row, sign * truss.web_perm[j] / limit))
        for index, chord in enumerate(truss.chords):
            co = {key: power(val.T) for key, val in states[index].items()}
            member = chord["member"]
            modes = [(n / (member["A"] * stress), m * member["c"] / (member["I"] * stress), 0.)
                     for n in (-1., 1.) for m in (-1., 1.)]
            if index == 0:
                modes.append((0., 0., -1. / truss.sag_limit))
            for panel in range(len(chord["lengths"])):
                for axial, bend, sag in modes:
                    self.labels.append((member, "SAG" if sag else "STRESS"))
                    self.fields.append(SectionField(truss, co, index, panel, axial, bend, sag, weight, stress, modulus))

    def evaluate(self, positions):
        """Attained upper capacity (kN) and witnessing configurations. The relative enclosure is RELATIVE_GAP."""
        s = np.asarray(positions, dtype=float)
        joints = self.truss.xb
        d, far = self.clearance, s + self.wheelbase
        locations = np.column_stack([s, far, s-d, s+d, far-d, far+d])
        panels = np.clip(np.searchsorted(joints, locations, side="right") - 1, 0, len(joints)-2)
        frac = (locations - joints[panels]) / self.truss.lb[panels]
        best = np.full(len(s), np.inf)
        detail = np.zeros((len(s), 4))
        all_ids = np.arange(len(s))
        span = joints[-1]
        weight, heavy = self.weight, self.heavy

        def sample(row, ids):
            return row_samples(row, ids, locations, panels, frac, d, span)

        for code, (row, dead) in enumerate(self.web_lines):
            endpoints(best, detail, all_ids, code, 0., dead, sample(row, all_ids), weight, heavy)

        def cover(field, code, a, b, ids, ra, rb, da, db, sa, sb):
            geometry = (field.beta, field.curvature, field.origin, field.length, field.diagonal_curvature)
            sub = inspect_cell(ids, best, detail, locations, joints, code, (a, b), (da, db), geometry,
                               field.diagonal, (ra, rb), (sa, sb), weight, d, heavy)
            if len(sub):
                active = ids[sub]
                mid = .5 * (a + b)
                if not a < mid < b:
                    raise ArithmeticError("section capacity enclosure did not converge")
                rm, dm = field.row(mid)
                sm = sample(rm, active)
                endpoints(best, detail, active, code, mid, dm, sm, weight, heavy)
                cover(field, code, a, mid, active, ra, rm, da, dm, sa[sub], sm)
                cover(field, code, mid, b, active, rm, rb, dm, db, sm, sb[sub])

        # Find attained endpoint bounds first. They make the later spatial search local.
        pending = []
        for code, field in enumerate(self.fields, len(self.web_lines)):
            a, b = field.origin, field.origin + field.length
            ra, da = field.row(a)
            rb, db = field.row(b)
            sa, sb = sample(ra, all_ids), sample(rb, all_ids)
            endpoints(best, detail, all_ids, code, a, da, sa, weight, heavy)
            endpoints(best, detail, all_ids, code, b, db, sb, weight, heavy)
            if field.curvature:
                pending.append((code, field, a, b, ra, rb, da, db, sa, sb))
        for code, field, a, b, ra, rb, da, db, sa, sb in pending:
            cover(field, code, a, b, all_ids, ra, rb, da, db, sa, sb)
        return best, detail
