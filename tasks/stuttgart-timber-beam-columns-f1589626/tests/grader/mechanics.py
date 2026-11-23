"""Piecewise equilibrium of the pin-ended member with both sectional stress traces."""

import csv
from decimal import Decimal, localcontext
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np


LENGTH = 8100
STEP = 150
MOMENT = 120 * 480**2 // 6
INERTIA = 120 * 480**3 // 12
PRELOAD = 68000
RELEASE_RATIO = Decimal("1.0009")


def read_stock(directory):
    """Keep the published decimal material values available for precision refinement."""
    with (Path(directory) / "data_beam_0.txt").open() as file:
        rows = csv.reader(file)
        names = [name.strip() for name in next(rows)]
        records = [dict(zip(names, row)) for row in rows if row]
    joints = np.array([int(r["finger joint"]) != 0 for r in records])

    def strength(name):
        return np.array([r[name + "_FJ" if joint else name].strip()
                         for r, joint in zip(records, joints)])

    return {
        "segment": np.array([int(r["n segment"]) for r in records]),
        "joint": joints,
        "ft": strength("ft"),
        "fc": strength("fc"),
        "modulus": np.array([[r["Et_le"].strip(), r["Et_ri"].strip()] for r in records]),
    }


def read_survey(directory):
    """Keep the original decimals, with full and half stations in physical order."""
    with (Path(directory) / "data_Imp_meas_0.txt").open() as file:
        rows = [np.array(r) for r in csv.reader(file) if r]

    def interleave(full, half):
        values = np.empty(full.size + half.size, dtype=object)
        values[::2] = full
        values[1::2] = half
        return values

    return interleave(*rows[:2]), interleave(*rows[2:])


def interp_stations(values, half_mm, number):
    """Interpolate at exact integer half-mm positions on the 40.5-mm survey grid."""
    index = np.minimum(half_mm // 81, len(values) - 2)
    weight = (half_mm - 81 * index) / number(81)
    return values[index] + weight * (values[index + 1] - values[index])


def inverse2(a):
    b = a.copy()
    b[..., 0, 0], b[..., 1, 1] = a[..., 1, 1], a[..., 0, 0]
    b[..., 0, 1], b[..., 1, 0] = -a[..., 0, 1], -a[..., 1, 0]
    det = a[..., 0, 0] * a[..., 1, 1] - a[..., 0, 1] * a[..., 1, 0]
    return b / np.asarray(det)[..., None, None]


def sections(modulus, angles, number):
    """Integrate the force/strain matrix using exact rectangle moments in material axes."""
    ys, zs = np.meshgrid([-30, 30], np.arange(-220, 240, 40), indexing="ij")
    y, z = ys.ravel(), zs.ravel()
    area_modulus = (60 * 40) * modulus
    a = area_modulus.sum(axis=1)
    b, c = area_modulus @ y, area_modulus @ z
    d = area_modulus @ (y**2 + 60**2 / number(12))
    e = area_modulus @ (y * z)
    f = area_modulus @ (z**2 + 40**2 / number(12))
    # The matrix for generalized section resultants in local coordinates is
    # [[a,b,c], [b,d,e], [c,e,f]]. Invert by its six distinct cofactors.
    aa, ab, ac = d*f-e*e, c*e-b*f, b*e-c*d
    bb, bc, cc = a*f-c*c, b*c-a*e, a*d-b*b
    determinant = a*aa + b*ab + c*ac
    local = np.stack((aa, ab, ac, ab, bb, bc, ac, bc, cc), axis=-1).reshape(-1, 3, 3)
    local /= determinant[:, None, None]
    if number is np.longdouble:
        c, s = np.cos(angles), np.sin(angles)
    else:
        c = np.array([mp.cos(v) for v in angles])
        s = np.array([mp.sin(v) for v in angles])
    transform = np.zeros_like(local)
    transform[:, 0, 0] = 1
    transform[:, 1:, 1:] = np.stack((c, -s, s, c), axis=-1).reshape(-1, 2, 2)
    # [e, -u_y'', -u_z''] = compliance [N, M_y, M_z].
    strain_from_load = local @ transform.transpose(0, 2, 1)
    curvature = (transform @ strain_from_load)[:, 1:, :]
    offsets = np.array([[-30, -20], [-30, 20], [30, -20], [30, 20]])
    corners = (np.stack((y, z), axis=-1)[:, None, :] + offsets).reshape(96, 2)
    basis = np.column_stack((np.ones(len(corners), dtype=int), corners))
    stress = (basis @ strain_from_load) * np.repeat(modulus, 4, axis=1)[..., None]
    return curvature, stress


def interval_loads(curvature, data, number, degree):
    """Use exact element end relations and impose slope equilibrium at every breakpoint."""
    # These coordinates in half-mm units include every change in the ODE coefficients.
    x = np.union1d(np.arange(55) * (2 * STEP), np.arange(len(data)) * 81)
    h = (np.diff(x) / number(2))[:, None, None]
    section = ((x[1:] + x[:-1]) // (4 * STEP)).astype(int)
    bow = interp_stations(data, x, number)
    count = len(section)
    dtype = curvature.dtype
    loads = np.zeros((count, 3, 2), dtype=dtype)
    loads[:, 0, 0] = -PRELOAD
    loads[:, 1, 0] = PRELOAD * bow[:-1]
    loads[:, 2, 1] = MOMENT
    ramp = np.zeros_like(loads)
    ramp[:, 1, 0] = PRELOAD * np.diff(bow)
    k = -h**2 * PRELOAD * curvature[section, :, 1:]

    # In local t=(x-x0)/h the solution is C(t) u_left + S(t) h u'_left + F(t).
    # Each of the two load columns has affine forcing. The matrix power series
    # solves the coupled equation directly, without diagonalizing section tensors.
    eye = np.eye(2, dtype=dtype)
    coef = np.zeros((degree + 1, count, 2, 6), dtype=dtype)
    coef[0, :, :, :2], coef[1, :, :, 2:4] = eye, eye
    for j in range(2, degree + 1):
        rhs = k @ coef[j - 2]
        if j in (2, 3):
            rhs[:, :, 4:] -= h**2 * (curvature[section] @ (loads if j == 2 else ramp))
        coef[j] = rhs / number(j * (j - 1))
    end = coef.sum(axis=0)
    derivative = (np.arange(1, degree + 1)[:, None, None, None] * coef[1:]).sum(axis=0)
    c, s, force = end[:, :, :2], end[:, :, 2:4], end[:, :, 4:]
    dc, ds, df = derivative[:, :, :2], derivative[:, :, 2:4], derivative[:, :, 4:]
    sinv = inverse2(s)
    ll = sinv @ c / h
    lr = -sinv / h
    rr = ds @ sinv / h
    # Element end-slope vector is [-u'_left,+u'_right]. Continuity of
    # physical slopes produces a two-by-two block tridiagonal Dirichlet system.
    rl = (dc - ds @ sinv @ c) / h
    fleft = sinv @ force / h
    fright = (df - ds @ sinv @ force) / h
    diag = rr[:-1] + ll[1:]
    rhs = -fright[:-1] - fleft[1:]
    for j in range(1, count - 1):
        factor = rl[j] @ inverse2(diag[j - 1])
        diag[j] -= factor @ lr[j]
        rhs[j] -= factor @ rhs[j - 1]
    displacement = np.zeros((count + 1, 2, 2), dtype=dtype)
    for j in reversed(range(count - 1)):
        displacement[j + 1] = inverse2(diag[j]) @ (rhs[j] - lr[j + 1] @ displacement[j + 2])

    slope = sinv @ (displacement[1:] - c @ displacement[:-1] - force)
    solution = coef[..., :2] @ displacement[:-1] + coef[..., 2:4] @ slope + coef[..., 4:]
    load_poly = np.zeros((degree + 1, count, 3, 2), dtype=dtype)
    load_poly[:, :, 1:, :] = PRELOAD * solution
    load_poly[0] += loads
    load_poly[1] += ramp
    return section, load_poly


def bernstein_coefficients(power, number):
    """Convert a matrix-valued Taylor polynomial on [0,1] to Bernstein coefficients."""
    degree = len(power) - 1
    coeff = []
    for i in range(degree + 1):
        value = power[0].copy()
        for j in range(1, i + 1):
            value += power[j] * (number(comb(i, j)) / comb(degree, j))
        coeff.append(value)
    return np.stack(coeff, axis=1)


def divide_bounds(stress, ft, fc):
    preload, bending = stress[..., 0], stress[..., 1]
    reserve = np.minimum(ft[:, None] - preload, fc[:, None] + preload)
    ratio = np.full_like(preload, np.inf)
    np.divide(reserve, abs(bending), out=ratio, where=bending != 0)
    # Any negative preload reserve gives a zero bound, even at zero bending sensitivity.
    return np.maximum(0, np.where(reserve < 0, 0, ratio))


def split_curve(coeff):
    """de Casteljau subdivision at 1/2 preserves the entire curve in both children."""
    work = coeff.copy()
    left, right = np.empty_like(work), np.empty_like(work)
    degree = coeff.shape[1] - 1
    left[:, 0], right[:, degree] = work[:, 0], work[:, -1]
    for step in range(1, degree + 1):
        work = (work[:, :-1] + work[:, 1:]) * 0.5
        left[:, step], right[:, degree - step] = work[:, 0], work[:, -1]
    return np.concatenate((left, right))


def envelope(stress, ft, fc, relative_gap):
    # At either sign of moment, a nonnegative strength margin at every control
    # point guarantees the same throughout the interval (convex-hull property).
    ratios = divide_bounds(stress, ft, fc)
    upper = ratios[:, [0, -1]].min()
    lower = ratios.min(axis=1)
    active = np.flatnonzero(lower <= upper)
    while len(active):
        floor = lower[active].min()
        if upper - floor <= relative_gap * max(1, upper):
            break
        stress = split_curve(stress[active])
        ft, fc = np.tile(ft[active], 2), np.tile(fc[active], 2)
        ratios = divide_bounds(stress, ft, fc)
        upper = min(upper, ratios[:, [0, -1]].min())
        lower = ratios.min(axis=1)
        active = np.flatnonzero(lower < upper)
    return upper


def material(data, number):
    if number is np.longdouble:
        return np.asarray(data, dtype=number)
    return np.vectorize(lambda value: number(str(value)), otypes=[object])(data)


def calculate(layout, stock, survey, number, degree, precision, possible=None):
    row, side = layout[..., 0], layout[..., 1]
    modulus = material(stock["modulus"][row, side].reshape(54, 24), number)
    bow, angles = (material(values, number) for values in survey)
    twist = interp_stations(angles, (2 * np.arange(54) + 1) * STEP, number)
    curvature, stress_map = sections(modulus, twist, number)
    apparent = number(LENGTH) / (INERTIA * STEP * curvature[:, 1, 2].sum())
    section, loads = interval_loads(curvature, 1000 * bow, number, degree)
    b = bernstein_coefficients(loads, number)
    # Each width/layer rectangle contributes all four one-sided corner traces.
    if possible is None:
        traces = np.arange(len(section) * 96)
    else:
        traces = possible
    interval, corner = np.divmod(traces, 96)
    record = row.reshape(54, 24)[section[interval], corner // 4]
    ft = material(stock["ft"][record], number)
    fc = material(stock["fc"][record], number)
    if possible is None:
        stress = (stress_map[section, None, :, :] @ b).transpose(0, 2, 1, 3)
        stress = stress.reshape(len(traces), degree + 1, 2)
    else:
        stress = (stress_map[section[interval], corner][:, None, None, :] @ b[interval])[:, :, 0]
    capacity = envelope(stress, ft, fc, precision)
    # Keep all potential winning traces within a forward-roundoff guard when
    # reducing the expensive arbitrary-precision calculation. It does not alter
    # a tolerance or a pass threshold. The band is deliberately conservative
    # for the long-double block solve with these fixed stocks and pinned spans.
    if possible is None:
        lower = divide_bounds(stress, ft, fc).min(axis=1)
        possible = np.flatnonzero(lower <= capacity + number("1e-7"))
    return (capacity, apparent), possible


class MemberRating:
    """Take each rating's minimum over installations at the comparison precision."""

    def __init__(self, layout, stock, survey):
        # The turn maps fixture station and width indices to their opposites.
        # Depth and the (stock row, source half) carried by a cell stay intact.
        # Both installed maps then see the same fixture survey and twist sign.
        self.installations = (layout, layout[::-1, ::-1, :, :])
        self.stock, self.survey = stock, survey
        self.values, self.possible = self.ratings(np.longdouble, 14,
                                                 8 * np.finfo(np.longdouble).eps,
                                                 (None, None))
        self.refined = {}

    def ratings(self, number, degree, precision, possible):
        results = [calculate(installed, self.stock, self.survey, number, degree, precision, traces)
                   for installed, traces in zip(self.installations, possible)]
        (forward, forward_traces), (turned, turned_traces) = results
        minima = (min(forward[0], turned[0]), min(forward[1], turned[1]))
        return minima, (forward_traces, turned_traces)

    def __getitem__(self, index):
        return self.values[index]

    def compare(self, field, value, tolerance=None):
        supplied = Decimal(str(value))

        def verdict(estimate, guard):
            low, high = estimate - guard, estimate + guard
            if tolerance is None:
                if high < supplied:
                    return False
                if low >= supplied:
                    return True
            else:
                if supplied < low - tolerance or supplied > high + tolerance:
                    return False
                if high - tolerance <= supplied <= low + tolerance:
                    return True
            return None

        # Use a wider interval to trigger recomputation if roundoff might
        # cross a decision boundary. It never expands the acceptance band.
        estimate = Decimal(str(self.values[field]))
        decision = verdict(estimate, Decimal("1e-8") * max(1, abs(estimate)))
        if decision is not None:
            return decision
        digits = 42
        while True:
            with mp.workdps(digits), localcontext() as ctx:
                ctx.prec = digits
                if digits not in self.refined:
                    # On the fixed spans the small matrix Taylor coefficients
                    # fall factorially. Increase both order and precision.
                    degree = max(22, digits // 2 + 2)
                    values, _ = self.ratings(mp.mpf, degree, mp.power(10, 7 - digits),
                                             self.possible)
                    self.refined[digits] = values
                estimate = Decimal(mp.nstr(self.refined[digits][field], digits))
                guard = Decimal(1).scaleb(10 - digits) * max(1, abs(estimate))
                decision = verdict(estimate, guard)
                if decision is not None:
                    return decision
                if estimate == 0 and field == 0:
                    # A negative preload margin was detected before bending.
                    return verdict(estimate, Decimal(0))
            digits *= 2


def rate_member(layout, stock, survey):
    """Rate the nominal stock under the proof load and pure bending in both installations."""
    return MemberRating(layout, stock, survey)


def within_report(rating, field, reported):
    """Apply the stated absolute reporting tolerance to the original number."""
    tol = (Decimal("0.002"), Decimal("0.2"))[field]
    return rating.compare(field, reported, tol)


def meets_target(rating, field, target):
    """Compare the recomputed rating to the grade target with its release margin."""
    return rating.compare(field, RELEASE_RATIO * Decimal(str(target)))
