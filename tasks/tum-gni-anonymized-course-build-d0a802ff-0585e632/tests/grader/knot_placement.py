"""Places the knots of a piecewise-linear curve inside a corridor.

The corridor is given at stations sig[0] < ... < sig[m] by a lower and an upper
ordinate and is straight in between.  The curve has knots (x[i], y[i]) with
increasing x.  Every knot may be moved along the station axis by up to `reach`,
each by its own amount.  The question answered here is whether the knots can be
placed at stations t[i] such that

  * |t[i] - x[i]| <= reach for every knot,
  * t[0] <= t[1] <= ... (the knots keep their order),
  * t[0] <= sig[0] and t[-1] >= sig[m] (the curve covers the corridor), and
  * the curve that is linear between the points (t[i], y[i]) stays inside the
    corridor at every station of sig[0]..sig[m].

Nothing is sampled.  With t and t' the stations of two consecutive knots of
ordinates a and b, the ordinate of the curve at a corridor station s between
them is (a (t' - s) + b (s - t)) / (t' - t), so a corridor bound at s is a
half-plane in (t, t').  Because curve and corridor are both straight between
corridor stations, the curve is inside the corridor between two knots exactly
when it is inside at the corridor stations between them and at the knots
themselves.  As long as t and t' stay in fixed corridor cells the set of
admissible (t, t') is therefore a convex polygon, and the set of stations that
a knot can take, given everything before it, is a union of intervals.  These
sets are carried from the first knot to the last one.

Corridor cell c is sig[c]..sig[c + 1]; cell -1 is everything up to sig[0] and
cell m everything from sig[m] on.  Outside sig[0]..sig[m] there is no corridor.
"""
import numpy as np

EDGE = 1e-12        # round-off allowed on the edge of a half-plane, in units of the station axis
THIN = 1e-9         # a set of stations may be this much narrower than nothing through round-off


def _clip(px, py, n, ca, cb, cc, qx, qy):
    """Part of a convex polygon where ca * px + cb * py + cc >= 0, written to (qx, qy); returns its vertex count.

    The result has at most 2 n vertices, which (qx, qy) must be able to hold."""
    tol = EDGE * (abs(ca) + abs(cb))
    m = 0
    kp = n - 1
    fp = ca * px[kp] + cb * py[kp] + cc
    for k in range(n):
        f = ca * px[k] + cb * py[k] + cc
        if (f >= -tol) != (fp >= -tol):               # the edge from vertex kp to vertex k crosses the line
            r = fp / (fp - f)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            vx = px[kp] + r * (px[k] - px[kp])
            vy = py[kp] + r * (py[k] - py[kp])
            if m == 0 or qx[m - 1] != vx or qy[m - 1] != vy:
                qx[m] = vx
                qy[m] = vy
                m += 1
        if f >= -tol and (m == 0 or qx[m - 1] != px[k] or qy[m - 1] != py[k]):
            qx[m] = px[k]
            qy[m] = py[k]
            m += 1
        kp = k
        fp = f
    return m


def _cell_stretch(cp, b, xb, sig, low, upp, reach):
    """Stations of corridor cell cp where a knot of ordinate b is inside the corridor and within reach of xb.

    Returned as offsets from xb; the first is larger than the second when there are none."""
    m = sig.size - 1
    if cp < 0:
        return -reach, min(reach, sig[0] - xb)
    if cp >= m:
        return max(-reach, sig[m] - xb), reach
    first, last = 0.0, 1.0                     # as a share of the cell
    rise = low[cp + 1] - low[cp]
    if rise > 0.0:
        last = min(last, (b - low[cp]) / rise)
    elif rise < 0.0:
        first = max(first, (low[cp] - b) / (-rise))
    elif b < low[cp]:
        return 1.0, -1.0
    rise = upp[cp + 1] - upp[cp]
    if rise > 0.0:
        first = max(first, (b - upp[cp]) / rise)
    elif rise < 0.0:
        last = min(last, (upp[cp] - b) / (-rise))
    elif b > upp[cp]:
        return 1.0, -1.0
    if first > last:
        return 1.0, -1.0
    width = sig[cp + 1] - sig[cp]
    return max(-reach, sig[cp] + first * width - xb), min(reach, sig[cp] + last * width - xb)


def _sweep(x, y, sig, low, upp, reach):
    """Stations every knot can take, given the knots before it.

    Returns (stuck, first, lo, hi, cell).  The sets of knot i are the intervals
    first[i]..first[i + 1] - 1 of lo, hi (offsets from x[i]) and cell (their corridor
    cell), ordered by cell.  `stuck` is the first knot without any station, or -1."""
    n = x.size
    m = sig.size - 1
    first = np.zeros(n + 1, np.int64)
    size = 8 * n + 64
    lo = np.empty(size)
    hi = np.empty(size)
    cell = np.empty(size, np.int64)
    count = 0
    room = 64                                   # vertices that the polygon buffers hold
    px = np.empty(room)
    py = np.empty(room)
    qx = np.empty(room)
    qy = np.empty(room)
    wl = np.empty(64)
    wh = np.empty(64)

    # the first knot stands at or before sig[0]
    al = -reach
    be = min(reach, sig[0] - x[0])
    if al <= be:
        lo[0] = al
        hi[0] = be
        cell[0] = -1
        count = 1
    first[1] = count
    if count == 0:
        return 0, first, lo, hi, cell

    for i in range(n - 1):
        a = y[i]
        b = y[i + 1]
        xa = x[i]
        xb = x[i + 1]
        gap = xb - xa
        s0 = first[i]
        s1 = first[i + 1]
        if wl.size < s1 - s0:
            wl = np.empty(2 * (s1 - s0))
            wh = np.empty(2 * (s1 - s0))
        c1 = np.searchsorted(sig, xb - reach, side="left") - 1      # first cell that reaches the knot
        c2 = np.searchsorted(sig, xb + reach, side="right") - 1     # last cell that starts within its reach
        if c1 < cell[s0]:
            c1 = cell[s0]                                           # not before the knot behind
        if i + 2 == n and c1 < m:
            c1 = m                                                  # the last knot stands at or beyond sig[m]
        for cp in range(c1, c2 + 1):
            al, be = _cell_stretch(cp, b, xb, sig, low, upp, reach)
            if al > be:
                continue
            found = 0
            for k in range(s0, s1):
                c = cell[k]
                if c > cp:
                    break
                px[0] = lo[k]
                py[0] = al
                px[1] = hi[k]
                py[1] = al
                px[2] = hi[k]
                py[2] = be
                px[3] = lo[k]
                py[3] = be
                cnt = 4
                if c == cp:                                         # same cell: keep the order of the two knots
                    cnt = _clip(px, py, cnt, -1.0, 1.0, gap, qx, qy)
                    px, qx = qx, px
                    py, qy = qy, py
                for j in range(c + 1, cp + 1):                      # corridor stations between the two knots
                    if cnt == 0:
                        break
                    if 4 * cnt > room:                              # two clips can follow, each may double
                        room = 8 * cnt
                        grown = np.empty(room)
                        grown[:cnt] = px[:cnt]
                        px = grown
                        grown = np.empty(room)
                        grown[:cnt] = py[:cnt]
                        py = grown
                        qx = np.empty(room)
                        qy = np.empty(room)
                    sj = sig[j]
                    bound = low[j]
                    if a < bound or b < bound:
                        cnt = _clip(px, py, cnt, bound - b, a - bound,
                                    (a - bound) * (xb - sj) + (b - bound) * (sj - xa), qx, qy)
                        px, qx = qx, px
                        py, qy = qy, py
                        if cnt == 0:
                            break
                    bound = upp[j]
                    if a > bound or b > bound:
                        cnt = _clip(px, py, cnt, b - bound, bound - a,
                                    (bound - a) * (xb - sj) + (bound - b) * (sj - xa), qx, qy)
                        px, qx = qx, px
                        py, qy = qy, py
                if cnt == 0:
                    continue
                t0 = py[0]
                t1 = py[0]
                for v in range(1, cnt):
                    if py[v] < t0:
                        t0 = py[v]
                    if py[v] > t1:
                        t1 = py[v]
                t0 = max(t0, al)
                t1 = min(t1, be)
                if t0 > t1:
                    continue
                # keep the intervals of this cell ordered by their start
                pos = found
                while pos > 0 and wl[pos - 1] > t0:
                    wl[pos] = wl[pos - 1]
                    wh[pos] = wh[pos - 1]
                    pos -= 1
                wl[pos] = t0
                wh[pos] = t1
                found += 1
            if found == 0:
                continue
            if count + found > size:
                size = 2 * size + found
                grown = np.empty(size)
                grown[:count] = lo[:count]
                lo = grown
                grown = np.empty(size)
                grown[:count] = hi[:count]
                hi = grown
                larger = np.empty(size, np.int64)
                larger[:count] = cell[:count]
                cell = larger
            lo[count] = wl[0]
            hi[count] = wh[0]
            cell[count] = cp
            for v in range(1, found):
                if wl[v] <= hi[count] + EDGE:
                    if wh[v] > hi[count]:
                        hi[count] = wh[v]
                else:
                    count += 1
                    lo[count] = wl[v]
                    hi[count] = wh[v]
                    cell[count] = cp
            count += 1
        first[i + 2] = count
        if count == s1:
            return i + 1, first, lo, hi, cell
    return -1, first, lo, hi, cell


def _trace(x, y, sig, low, upp, first, lo, hi, cell):
    """One placement out of the sets of _sweep, taken from the last knot backwards.

    Every knot is put in the middle of the stations that suit the knot after it.
    Returns the stations and the knot at which none was found, or -1."""
    n = x.size
    t = np.empty(n)
    best = -1
    width = -1.0
    for k in range(first[n - 1], first[n]):
        if hi[k] - lo[k] > width:
            width = hi[k] - lo[k]
            best = k
    if best < 0:
        return t, n - 1
    off = 0.5 * (lo[best] + hi[best])
    cp = cell[best]
    t[n - 1] = x[n - 1] + off
    for i in range(n - 2, -1, -1):
        a = y[i]
        b = y[i + 1]
        xa = x[i]
        xb = x[i + 1]
        tb = xb + off
        best = -1
        width = -THIN
        pick = 0.0
        for k in range(first[i], first[i + 1]):
            c = cell[k]
            if c > cp:
                break
            tl = lo[k]
            th = hi[k]
            if c == cp:
                th = min(th, tb - xa)
            for j in range(c + 1, cp + 1):
                if th < tl - 1.0:
                    break
                sj = sig[j]
                bound = low[j]
                if a < bound or b < bound:
                    ca = bound - b
                    rest = (a - bound) * (tb - sj) + (b - bound) * (sj - xa)
                    if ca > 0.0:
                        tl = max(tl, -rest / ca)
                    elif ca < 0.0:
                        th = min(th, rest / (-ca))
                    elif rest < 0.0:
                        th = tl - 2.0
                bound = upp[j]
                if a > bound or b > bound:
                    ca = b - bound
                    rest = (bound - a) * (tb - sj) + (bound - b) * (sj - xa)
                    if ca > 0.0:
                        tl = max(tl, -rest / ca)
                    elif ca < 0.0:
                        th = min(th, rest / (-ca))
                    elif rest < 0.0:
                        th = tl - 2.0
            if th - tl > width:
                width = th - tl
                best = k
                pick = 0.5 * (tl + th)
        if best < 0:
            return t, i
        off = pick
        cp = cell[best]
        t[i] = xa + off
    return t, -1


_PLAIN = {"_clip": _clip, "_cell_stretch": _cell_stretch, "_sweep": _sweep, "_trace": _trace}
_READY = []


def _prepare():
    """Compiles the four functions above on first use.

    They are written in the subset of Python that numba compiles and give the same
    results either way; without a usable compiler they simply run as they are."""
    if _READY:
        return
    try:
        from numba import njit
        for name, function in _PLAIN.items():
            globals()[name] = njit(cache=False)(function)
        # a first small case, so that a compiler that cannot do its work is noticed here
        probe = _sweep(np.array([-0.5, 0.5, 1.5]), np.ones(3), np.array([0.0, 1.0]), np.full(2, 0.9),
                       np.full(2, 1.1), 1.0)
        if probe[0] != -1 or _trace(np.array([-0.5, 0.5, 1.5]), np.ones(3), np.array([0.0, 1.0]),
                                    np.full(2, 0.9), np.full(2, 1.1), *probe[1:])[1] != -1:
            raise ArithmeticError("compiled sweep disagrees with itself")
        _READY.append("compiled")
    except Exception:
        globals().update(_PLAIN)
        _READY.append("plain")


def place(x, y, sig, low, upp, reach):
    """Places the knots (x, y) inside the corridor (sig, low, upp), each within `reach` of its station.

    Returns (stuck, t): `stuck` is -1 and t holds the stations of one placement, or
    `stuck` is the first knot that cannot be placed whatever the knots before it do
    and t is None."""
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    sig = np.ascontiguousarray(sig, dtype=np.float64)
    low = np.ascontiguousarray(low, dtype=np.float64)
    upp = np.ascontiguousarray(upp, dtype=np.float64)
    if x.size < 2 or sig.size < 2 or np.any(np.diff(x) <= 0.0) or np.any(np.diff(sig) <= 0.0):
        raise ValueError("stations must increase")
    _prepare()
    stuck, first, lo, hi, cell = _sweep(x, y, sig, low, upp, float(reach))
    if stuck >= 0:
        return int(stuck), None
    t, lost = _trace(x, y, sig, low, upp, first, lo, hi, cell)
    if lost >= 0:
        raise ArithmeticError("no placement was traced back from knot %d" % lost)
    return -1, t


def excess(t, x, y, sig, low, upp, reach):
    """How far a placement is from what `place` promises; nothing above round-off for a sound one.

    Returns the largest of: the distance of a knot beyond its reach, a decrease of the
    stations, the shortfall of the cover, and the distance of the curve outside the
    corridor as a share of the corridor ordinate (checked at every corridor station and
    at every knot inside the corridor, which is everywhere the two can part)."""
    t = np.asarray(t, dtype=float)
    n, m = t.size, sig.size - 1
    worst = max(float(np.max(np.abs(t - x)) - reach), float(np.max(-np.diff(t))),
                float(t[0] - sig[0]), float(sig[m] - t[-1]), 0.0)
    scale = max(float(np.max(upp)), 1e-300)

    def outside(value, lower, upper):
        return float(np.max(np.maximum(lower - value, value - upper) / np.maximum(upper, 1e-9 * scale)))

    # corridor stations: the knot at or behind each one and the next knot
    behind = np.searchsorted(t, sig, side="right") - 1
    behind[m] = np.searchsorted(t, sig[m], side="left") - 1     # at the far end the first knot there counts
    behind = np.clip(behind, 0, n - 2)
    span = t[behind + 1] - t[behind]
    share = np.clip(np.divide(sig - t[behind], span, out=np.zeros(m + 1), where=span > 0.0), 0.0, 1.0)
    worst = max(worst, outside(y[behind] + (y[behind + 1] - y[behind]) * share, low, upp))
    inner = (t > sig[0]) & (t < sig[m])
    if inner.any():
        worst = max(worst, outside(y[inner], np.interp(t[inner], sig, low), np.interp(t[inner], sig, upp)))
    return worst
