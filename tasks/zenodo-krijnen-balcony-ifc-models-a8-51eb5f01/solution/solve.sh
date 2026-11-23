#!/usr/bin/env bash
# Cast-in-place concrete take-off for the Tekla IFC export /app/inputs/a8.ifc.
# Writes the take-off program and its result for a8; nothing else goes under /app/out.
set -euo pipefail

OUT=/app/out
mkdir -p "$OUT"
find "$OUT" -mindepth 1 -delete          # stale copies

cat > "$OUT/pour_takeoff.py" <<'POUR_TAKEOFF_PY'
#!/opt/conda/bin/python
"""Cast-in-place concrete take-off for Tekla Structures IFC exports.

usage: /opt/conda/bin/python pour_takeoff.py <model.ifc> <result.json>

Parts
    Every IfcElement that carries a body is a part and is aggregated by one
    IfcElementAssembly.  It has the first class that applies:
    void (Name SPARING), annotation (material MISCELLANEOUS/DUMMY, ignored),
    poured (`Cast unit type` of the assembly's property set `Tekla Assembly`
    is Insitu and the part is not an IfcPile; its kind is its Name),
    other construction (the rest, piles included).  A body outside every
    assembly and a poured part whose Name is not one of the nine kinds are
    outside the published conventions and stop the run.

Geometry
    IfcOpenShell evaluates every part in its own placement frame (metres,
    whatever the unit of the file); the placement is applied here in double
    precision relative to a local origin, so a site 1000 km from the origin,
    any rotation and either length unit give the same solids, in a local
    millimetre frame.  Curved faces are faceted with a 0.2 mm chord.  Each
    part becomes a closed Manifold solid; a body that cannot be evaluated or
    closed stops the run.

Volumes
    body(k) = U poured(k) - U poured(earlier kinds) - U other - U voids,
    integrated over the height part by part: the horizontal section of a
    part less the sections of everything that ranks before it (voids, other
    construction, earlier kinds, earlier parts of its kind), by exact 2-D
    polygon booleans.  The height is cut at the levels of horizontal faces
    and at the corners of the tilted faces that reach into the box of the
    part: in between the section is constant where no tilted face is met
    and is integrated with Gauss panels where one is.

Surface
    The boundary of the poured body lies on faces of the parts, so every
    planar face of every part near poured concrete is a candidate, looked at
    from both sides.  The model is read at 1 mm resolution: what counts on
    either side of a face is the material PROBE mm away from it along its
    normal, i.e. beyond every gap, film, sheet or stack of them thinner than
    that.  Two planes parallel to the face, PROBE mm behind and PROBE mm in
    front of it, are cut through all neighbouring parts and the material on
    each plane (kind k / other construction / nothing) is found with exact
    2-D polygon booleans: voids erase everything, other construction
    displaces concrete, the earliest kind of the pour sequence owns
    overlapping concrete.
    A point of a face is boundary of kind k when kind k lies behind it and
    no poured concrete lies in front of it (concrete in front makes an
    internal joint); it is `contact` when other construction lies in front,
    otherwise free: top / soffit / side by the outward normal.
    All faces of one stack of thin layers see the same two materials, so a
    stack is counted once, on a face at which the concrete ends (concrete
    FINE mm behind it, none FINE mm in front of it): the first such face
    among those lying within PROBE mm of one another.  A face with concrete
    on both sides at FINE mm is a cut through the concrete and never counts.
    Boundary that no such face stands for (stacks of layers thinner than
    FINE) is counted afterwards from the remaining faces, again once per
    stack.
    PROBE sits just under the 1 mm resolution, in the band the exports leave
    empty (no layer between 0.8 mm and 5 mm), so thinner layers can neither
    create nor hide surface, wherever the model sits and however it is
    rotated.
"""
import sys
import json
import time

import numpy as np
import ifcopenshell
import ifcopenshell.geom
import manifold3d as m3
from manifold3d import Manifold, CrossSection, FillRule, OpType
from rtree import index as rtree_index

KINDS = ['POER', 'FUNDERINGSBALK', 'VLOER', 'OPSTORTING', 'KOLOM', 'KOLOMKOP',
         'WAND', 'BALK', 'DRUKLAAG']
VOID_NAME = 'SPARING'
DUMMY_MATERIAL = 'MISCELLANEOUS/DUMMY'
PSET_ASSEMBLY = 'Tekla Assembly'
PROP_CAST = 'Cast unit type'
CAST_POURED = 'Insitu'

PROBE = 0.99         # mm: resolution; distance of the two probing planes from a face
FINE = 0.01          # mm: distance at which a face is tested for concrete on both sides
DEFLECTION = 0.0002  # m : chordal error used to facet curved faces (0.2 mm)
ANGULAR = 0.2        # rad: angular deflection of the mesher
FINE_PRECISION = 1e-7  # m : kernel precision for bodies it rejects at its default (0.01 mm)
PLANE_TOL = 2e-3     # mm: flatness tolerance when grouping triangles into faces
STEEP = 1.75e-4      # |nz| of a face 0.01 deg off the vertical: a face beyond that is inclined

SOLID_TYPES = ('SweptSolid', 'Brep', 'CSG', 'Clipping', 'AdvancedSweptSolid', 'AdvancedBrep',
               'SolidModel', 'SurfaceModel', 'Tessellation', 'MappedRepresentation')

POURED, VOID, OTHER, ANNOTATION = 'poured', 'void', 'other', 'annotation'
TRACE = None         # optional list collecting every counted patch (diagnostics)


class ModelError(Exception):
    """the model breaks a convention the take-off relies on."""


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------
def _cast_unit_type(assembly):
    for rel in getattr(assembly, 'IsDefinedBy', None) or ():
        if not rel.is_a('IfcRelDefinesByProperties'):
            continue
        pdef = rel.RelatingPropertyDefinition
        if pdef is None or not pdef.is_a('IfcPropertySet') or pdef.Name != PSET_ASSEMBLY:
            continue
        for prop in pdef.HasProperties or ():
            if prop.Name == PROP_CAST and prop.is_a('IfcPropertySingleValue'):
                val = prop.NominalValue
                if val is not None:
                    return str(val.wrappedValue)
    return None


def _material_names(mat, out, depth=0):
    if mat is None or depth > 6:
        return
    if mat.is_a('IfcMaterial'):
        out.add(mat.Name)
        return
    for attr in ('ForLayerSet', 'ForProfileSet', 'MaterialLayers', 'MaterialProfiles',
                 'MaterialConstituents', 'Materials', 'Material'):
        val = getattr(mat, attr, None)
        if val is None:
            continue
        for item in (val if isinstance(val, (tuple, list)) else (val,)):
            _material_names(item, out, depth + 1)


def _materials(element):
    names = set()
    for rel in getattr(element, 'HasAssociations', None) or ():
        if rel.is_a('IfcRelAssociatesMaterial'):
            _material_names(rel.RelatingMaterial, names)
    return names


def _has_body(element):
    rep = element.Representation
    for r in (rep.Representations or ()) if rep is not None else ():
        if r.RepresentationIdentifier == 'Body' or (
                r.RepresentationIdentifier is None and r.RepresentationType in SOLID_TYPES):
            return True
    return False


def classify(model):
    """-> {element id: (class, kind index, name)} for every element with a body."""
    assemblies = {}
    for rel in model.by_type('IfcRelAggregates'):
        whole = rel.RelatingObject
        if whole is None or not whole.is_a('IfcElementAssembly'):
            continue
        for child in rel.RelatedObjects or ():
            assemblies.setdefault(child.id(), {})[whole.id()] = whole
    cast_cache = {}

    def cast(assembly):
        key = assembly.id()
        if key not in cast_cache:
            cast_cache[key] = _cast_unit_type(assembly)
        return cast_cache[key]

    parts = {}
    for el in model.by_type('IfcElement'):
        if el.is_a('IfcFeatureElement') or not _has_body(el):
            continue                              # openings are applied by the kernel; no body: no part
        name = el.Name or ''
        owners = list(assemblies.get(el.id(), {}).values())
        if len(owners) != 1:
            raise ModelError('#%d %s %r is aggregated by %d assemblies, expected one'
                             % (el.id(), el.is_a(), name, len(owners)))
        if name == VOID_NAME:
            parts[el.id()] = (VOID, -1, name)
        elif _materials(el) == {DUMMY_MATERIAL}:
            parts[el.id()] = (ANNOTATION, -1, name)
        elif cast(owners[0]) == CAST_POURED and not el.is_a('IfcPile'):
            if name not in KINDS:
                raise ModelError('#%d %s: poured part of unknown kind %r'
                                 % (el.id(), el.is_a(), name))
            parts[el.id()] = (POURED, KINDS.index(name), name)
        else:
            parts[el.id()] = (OTHER, -1, name)
    return parts


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
class Part:
    __slots__ = ('eid', 'cls', 'kind', 'name', 'man', 'lo', 'hi')


def _weld(verts, tris, tol):
    key = np.round(verts / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.reshape(-1)
    t = inv[tris]
    ok = (t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])
    return verts[first], t[ok]


def _manifold(verts, tris):
    if len(tris) < 4:
        return None
    mesh = m3.Mesh64(vert_properties=np.ascontiguousarray(verts, dtype=np.float64),
                     tri_verts=np.ascontiguousarray(tris, dtype=np.uint64))
    man = Manifold(mesh)
    if man.status() != m3.Error.NoError or man.is_empty():
        return None
    return man


def _shells(tris, nv):
    """split a triangle list into vertex-connected shells."""
    par = np.arange(nv)

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for a, b, c in tris:
        ra = find(a)
        par[find(b)] = ra
        par[find(c)] = ra
    root = np.array([find(a) for a in tris[:, 0]])
    return [tris[root == r] for r in np.unique(root)]


def _resolve(man):
    """Manifold made of several shells -> proper solid (overlapping shells are
    united, inside-out shells are cavities)."""
    pieces = man.decompose()
    if len(pieces) <= 1:
        return man if man.volume() > 0 else None
    solid = [p for p in pieces if p.volume() > 0]
    holes = [p for p in pieces if p.volume() < 0]
    if not solid:
        return None
    out = Manifold.batch_boolean(solid, OpType.Add)
    if holes:
        flipped = []
        for h in holes:
            mesh = h.to_mesh64()
            hm = _manifold(np.asarray(mesh.vert_properties)[:, :3],
                           np.asarray(mesh.tri_verts)[:, ::-1])
            if hm is not None:
                flipped.append(hm)
        if flipped:
            out = out - Manifold.batch_boolean(flipped, OpType.Add)
    return out


def _close_t_junctions(v, t, tol):
    """split triangle edges at vertices lying on them: a shell whose needle
    faces were dropped by the kernel is closed in space but not edge to edge."""
    t = [tuple(int(i) for i in tri) for tri in t]
    for _ in range(200):
        count = {}
        for a, b, c in t:
            for e in ((a, b), (b, c), (c, a)):
                count[e] = count.get(e, 0) + 1
        loose = {e for e, k in count.items() if count.get((e[1], e[0]), 0) < k}
        if not loose:
            break
        ends = np.array(sorted({i for e in loose for i in e}))
        out, changed = [], False
        for tri in t:
            split = None
            for j in range(3):
                a, b, x = tri[j], tri[(j + 1) % 3], tri[(j + 2) % 3]
                if (a, b) not in loose:
                    continue
                ab = v[b] - v[a]
                length = float(np.linalg.norm(ab))
                if length <= 2.0 * tol:
                    continue
                w = v[ends] - v[a]
                s = (w @ ab) / length
                off = np.linalg.norm(w - np.outer(s / length, ab), axis=1)
                on = (s > tol) & (s < length - tol) & (off < tol)
                if on.any():
                    split = (a, b, x, int(ends[on][int(np.argmin(s[on]))]))
                    break
            if split is None:
                out.append(tri)
            else:
                a, b, x, c = split
                out += [(a, c, x), (c, b, x)]
                changed = True
        t = out
        if not changed:
            break
    return np.array(t, dtype=np.int64).reshape(-1, 3)


def build_solid(verts, tris):
    """triangles of a closed surface (mm) -> Manifold, or None."""
    for heal in (False, True):
        for tol in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
            v, t = _weld(verts, tris, tol)
            if heal:
                t = _close_t_junctions(v, t, max(tol, 1e-4))
            man = _manifold(v, t)
            if man is None:
                shells = [_manifold(v, s) for s in _shells(t, len(v))]
                if not shells or any(s is None for s in shells):
                    continue
                man = Manifold.compose(shells)
            if man.volume() < 0 and len(man.decompose()) <= 1:
                man = _manifold(v, t[:, ::-1])
                if man is None:
                    continue
            man = _resolve(man)
            if man is not None and not man.is_empty():
                return man
    return None


def load_parts(path):
    model = ifcopenshell.open(path)
    classes = classify(model)
    counts = {}
    for c in classes.values():
        counts[c[0]] = counts.get(c[0], 0) + 1
    log('parts:', counts)
    wanted = [model.by_id(i) for i, c in classes.items() if c[0] != ANNOTATION]
    if not wanted:
        return []
    settings = ifcopenshell.geom.settings()
    settings.set('use-world-coords', False)
    settings.set('mesher-linear-deflection', DEFLECTION)
    settings.set('mesher-angular-deflection', ANGULAR)
    raw = {}

    def take(shape):
        geom = shape.geometry
        v = np.array(geom.verts, dtype=np.float64).reshape(-1, 3)
        t = np.array(geom.faces, dtype=np.int64).reshape(-1, 3)
        if len(v) and len(t):
            m = np.array(shape.transformation.matrix, dtype=np.float64).reshape(4, 4).T
            raw.setdefault(shape.id, []).append((v, t, m[:3, :3].copy(), m[:3, 3].copy()))

    it = ifcopenshell.geom.iterator(settings, model, 1, include=wanted)
    if it.initialize():
        while True:
            take(it.get())
            if not it.next():
                break
    settings.set('precision', FINE_PRECISION)
    for el in wanted:                      # bodies thinner than the kernel's default precision
        if el.id() not in raw:
            try:
                take(ifcopenshell.geom.create_shape(settings, el))
            except Exception:
                pass
    missing = [i for i in classes if classes[i][0] != ANNOTATION and i not in raw]
    if missing:
        raise ModelError('no body geometry for parts %s' % missing[:20])
    # local frame: metres -> millimetres around the poured concrete, with the
    # large translation removed before it meets the small coordinates
    frame_ids = [i for i in raw if classes[i][0] == POURED] or list(raw)
    lo = np.min([(v @ r.T + p).min(0) for i in frame_ids for v, _, r, p in raw[i]], axis=0)
    hi = np.max([(v @ r.T + p).max(0) for i in frame_ids for v, _, r, p in raw[i]], axis=0)
    origin = np.round((lo + hi) / 2.0)
    parts = []
    for eid in sorted(raw):
        cls, kind, name = classes[eid]
        mans = []
        for v, t, r, p in raw[eid]:
            man = build_solid((v @ r.T + (p - origin)) * 1000.0, t)
            if man is None:
                raise ModelError('#%d %s: the body is not a closed solid' % (eid, name))
            mans.append(man)
        part = Part()
        part.eid, part.cls, part.kind, part.name = eid, cls, kind, name
        part.man = mans[0] if len(mans) == 1 else Manifold.batch_boolean(mans, OpType.Add)
        box = part.man.bounding_box()
        part.lo = np.array(box[:3])
        part.hi = np.array(box[3:])
        parts.append(part)
    return parts


# --------------------------------------------------------------------------
# volumes
# --------------------------------------------------------------------------
def boxes_index(los, his, pad):
    prop = rtree_index.Property()
    prop.dimension = 3
    if len(los) == 0:
        return rtree_index.Index(properties=prop)
    return rtree_index.Index(((i, tuple(los[i] - pad) + tuple(his[i] + pad), None)
                              for i in range(len(los))), properties=prop)


def cs_union(items):
    items = [c for c in items if not c.is_empty()]
    if not items:
        return CrossSection()
    if len(items) == 1:
        return items[0]
    return CrossSection.batch_boolean(items, OpType.Add)


def z_levels(man):
    """-> (levels of the horizontal faces of a solid, its tilted triangles).
    A tilted triangle is neither horizontal nor vertical; its row holds its
    box (xlo, ylo, zlo, xhi, yhi, zhi), the level of its middle corner and
    whether it is tilted enough to matter for the choice of panels."""
    mesh = man.to_mesh64()
    v = np.asarray(mesh.vert_properties, dtype=np.float64)[:, :3]
    p = v[np.asarray(mesh.tri_verts, dtype=np.int64)]
    cr = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    dbl = np.linalg.norm(cr, axis=1)
    ok = dbl > 1e-9
    p = p[ok]
    nz = np.abs(cr[ok, 2]) / dbl[ok]
    flat = nz >= 1.0 - 1e-12
    hard = np.unique(p[flat][:, :, 2])
    sel = ~flat & (nz > 1e-9)
    t = p[sel]
    tilt = np.hstack([t.min(axis=1), t.max(axis=1), np.median(t[:, :, 2], axis=1)[:, None],
                      (nz[sel] > STEEP)[:, None]]).reshape(-1, 8)
    return hard, tilt


GAUSS = ((0.5 - 0.5 * 0.6 ** 0.5, 5.0 / 18.0), (0.5, 8.0 / 18.0), (0.5 + 0.5 * 0.6 ** 0.5, 5.0 / 18.0))


def volumes(parts):
    """volume of every kind (mm3): each poured part is integrated over the
    height, slice by slice, less what the parts ranking before it take away
    (other construction, voids, earlier kinds, earlier parts of its kind).
    Sections are exact 2-D polygon booleans.  The height is cut at the levels
    of horizontal faces and at the corners of every tilted triangle that
    reaches into the box of the part.  In an interval without such a triangle
    the section is constant; otherwise it is integrated by 3-point Gauss
    panels (exact while no edge of one section sweeps over a corner of
    another one, hence several panels wherever a face is really inclined)."""
    vols = [0.0] * len(KINDS)
    order = [i for i in sorted(range(len(parts)), key=lambda i: (parts[i].kind, i))
             if parts[i].cls == POURED]
    if not order:
        return vols
    place = {i: n for n, i in enumerate(order)}
    index = boxes_index([p.lo for p in parts], [p.hi for p in parts], -1e-7)
    levels = {}
    for i in order:
        part = parts[i]
        lo, hi = part.lo, part.hi
        near = []
        for j in index.intersection(tuple(lo) + tuple(hi)):
            if j == i:
                continue
            if parts[j].cls != POURED or place[j] < place[i]:
                near.append(j)
        for j in near + [i]:
            if j not in levels:
                levels[j] = z_levels(parts[j].man)
        z0, z1 = float(lo[2]), float(hi[2])
        tilt = {}                                 # tilted triangles reaching into the box of the part
        for j in near + [i]:
            t = levels[j][1]
            if len(t):
                t = t[(t[:, 0] < hi[0]) & (t[:, 3] > lo[0]) & (t[:, 1] < hi[1]) & (t[:, 4] > lo[1])
                      & (t[:, 2] < z1) & (t[:, 5] > z0)]
                if len(t):
                    tilt[j] = t
        hard = np.concatenate([levels[j][0] for j in near + [i]] + [np.array([z0, z1])])
        soft = np.unique(np.concatenate([t[:, (2, 5, 6)].ravel() for t in tilt.values()] + [np.zeros(0)]))
        soft = soft[(soft > z0) & (soft < z1)]
        if len(soft) > 4000:                      # very finely faceted: uniform intervals instead
            soft = np.linspace(z0, z1, 4001)[1:-1]
        dense = len(soft) > 400                   # short intervals: one panel in each
        cuts = np.unique(np.concatenate([hard[(hard >= z0) & (hard <= z1)], soft]))
        keep = [cuts[0]]
        for c in cuts[1:]:
            if c - keep[-1] > 1e-7:
                keep.append(c)
        total = 0.0
        for a, b in zip(keep[:-1], keep[1:]):
            h = b - a
            mid = 0.5 * (a + b)

            def tilted(j):
                """0: part j is prismatic in here, 1: all but, 2: it has inclined faces."""
                t = tilt.get(j)
                if t is None:
                    return 0
                m = (t[:, 2] < b) & (t[:, 5] > a)
                if not m.any():
                    return 0
                return 2 if t[m, 7].any() else 1

            active = [j for j in near if parts[j].lo[2] < b and parts[j].hi[2] > a]
            shape = {j: tilted(j) for j in active}
            flat = [j for j in active if shape[j] == 0]          # constant section in here
            bent = [j for j in active if shape[j] > 0]
            own = tilted(i)
            taken = cs_union([parts[j].man.slice(mid) for j in flat])
            base = None
            if own == 0:
                base = part.man.slice(mid) - taken
                if base.is_empty():
                    continue
                if not bent:
                    total += h * base.area()
                    continue

            def area(z):
                sec = base if base is not None else part.man.slice(z) - taken
                if bent:
                    sec = sec - cs_union([parts[j].man.slice(z) for j in bent])
                return sec.area()

            val = [area(a + x * h) for x, _ in GAUSS]
            panels = int(np.ceil(h / 50.0))
            if own == 2 or any(shape[j] == 2 for j in bent):
                if not dense:                     # an inclined face sweeps through: several panels
                    panels = max(panels, min(8, int(np.ceil(h / 0.25))))
            elif abs(val[1] - 0.5 * (val[0] + val[2])) * h <= 1e3:
                panels = 1                        # faces a hair off the vertical: linear in here
            elif not dense:
                panels = max(panels, 4)
            if panels == 1:
                total += h * sum(w * y for (_, w), y in zip(GAUSS, val))
                continue
            h /= panels
            for m in range(panels):
                total += h * sum(w * area(a + (m + x) * h) for x, w in GAUSS)
        vols[part.kind] += total
    return vols


# --------------------------------------------------------------------------
# surface
# --------------------------------------------------------------------------
class Face:
    __slots__ = ('part', 'n', 'd', 'tris', 'lo', 'hi')


def planar_faces(part_index, man):
    """group the triangles of a solid into planar faces."""
    mesh = man.to_mesh64()
    v = np.asarray(mesh.vert_properties, dtype=np.float64)[:, :3]
    p = v[np.asarray(mesh.tri_verts, dtype=np.int64)]          # (m, 3, 3)
    cr = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    dbl = np.linalg.norm(cr, axis=1)
    ok = dbl > 1e-9
    p, cr, dbl = p[ok], cr[ok], dbl[ok]
    nrm = cr / dbl[:, None]
    free = np.ones(len(p), dtype=bool)
    faces = []
    for i in np.argsort(-dbl, kind='stable'):
        if not free[i]:
            continue
        n = nrm[i]
        d = float(n @ p[i, 0])
        idx = np.nonzero(free & (nrm @ n > 1.0 - 1e-7))[0]
        if len(idx) > 1:
            idx = idx[np.abs(p[idx] @ n - d).max(axis=1) < PLANE_TOL]
        free[idx] = False
        f = Face()
        f.part = part_index
        w = dbl[idx]
        nn = (nrm[idx] * w[:, None]).sum(0)
        f.n = nn / np.linalg.norm(nn)
        f.tris = p[idx]
        f.d = float(((f.tris @ f.n).mean(axis=1) * w).sum() / w.sum())
        pts = f.tris.reshape(-1, 3)
        f.lo, f.hi = pts.min(0), pts.max(0)
        faces.append(f)
    return faces


def frame(n):
    e = np.zeros(3)
    e[int(np.argmin(np.abs(n)))] = 1.0
    u = np.cross(e, n)
    u /= np.linalg.norm(u)
    return u, np.cross(n, u)


def surface(parts, probe=PROBE, fine=FINE):
    """-> array (kinds, 4): top, side, soffit, contact in mm2."""
    nk = len(KINDS)
    res = np.zeros((nk, 4))
    poured = [p for p in parts if p.cls == POURED]
    if not poured:
        return res
    reach = probe + 0.5
    part_idx = boxes_index([p.lo for p in parts], [p.hi for p in parts], 0.0)
    pour_idx = boxes_index([p.lo for p in poured], [p.hi for p in poured], reach)

    rank = {POURED: 0, VOID: 1, OTHER: 2}
    faces = []
    for i in sorted(range(len(parts)), key=lambda i: (rank[parts[i].cls], i)):
        for f in planar_faces(i, parts[i].man):
            if pour_idx.count(tuple(f.lo) + tuple(f.hi)) > 0:
                faces.append(f)
    log('candidate faces:', len(faces))
    face_idx = boxes_index([f.lo for f in faces], [f.hi for f in faces], 0.0)

    def stood_for(fi, sign, frm, table, before):
        """what faces of the same stack already stand for on face fi, looking
        the same way (sign * normal), as polygons in the frame of face fi."""
        o, uv, qlo, qhi, corners, diag, box = frm
        n = faces[fi].n
        done = []
        for gi in face_idx.intersection(box):
            if gi == fi or (before and gi > fi) or not table[gi]:
                continue
            g = faces[gi]
            c = float(g.n @ n)
            if abs(c) < 0.5:
                continue
            polys = table[gi].get(sign if c > 0 else -sign)
            if not polys:
                continue
            grad = g.n @ uv
            c0 = float(g.n @ o - g.d)
            s = c0 + corners @ grad                      # distance of this face from plane g
            if s.min() >= probe or s.max() <= -probe:
                continue
            g2 = [(poly - o) @ uv for poly in polys]
            if c < 0:
                g2 = [poly[::-1] for poly in g2]
            if np.abs(s).max() < probe:
                done.extend(g2)
            else:                                        # planes cross: only the strip where they meet
                gn = float(np.linalg.norm(grad))
                e = grad / gn
                mid = (qlo + qhi) / 2.0
                mid = mid - e * ((c0 + mid @ grad) / gn)
                w = probe / gn
                along = np.array([-e[1], e[0]]) * diag
                strip = np.array([mid - along - e * w, mid + along - e * w,
                                  mid + along + e * w, mid - along + e * w])
                g_cs = CrossSection(g2, FillRule.Positive) ^ CrossSection([strip], FillRule.NonZero)
                done.extend(g_cs.to_polygons())
        return done

    def book(k, col, e, other, fi, sign):
        total = e.area()
        touch = (e ^ other).area()
        res[k, 3] += touch
        res[k, col] += total - touch
        if TRACE is not None:
            f = faces[fi]
            TRACE.append((k, col, total - touch, touch, parts[f.part].eid,
                          parts[f.part].name, (f.lo + f.hi) / 2.0, sign * f.n))

    ends = [None] * len(faces)      # per face and side: where it ends the concrete, 3-D polygons
    late = []                       # boundary seen from a face that does not end the concrete
    for fi, f in enumerate(faces):
        n = f.n
        box = tuple(f.lo - reach) + tuple(f.hi + reach)
        near = list(part_idx.intersection(box))
        if not any(parts[j].cls == POURED for j in near):
            continue
        u, v = frame(n)
        uv = np.array([u, v]).T
        o = (f.lo + f.hi) / 2.0
        o = o + n * (f.d - n @ o)                       # origin on the face plane
        f2 = (f.tris - o) @ uv
        patch = CrossSection(list(f2), FillRule.Positive)
        q = f2.reshape(-1, 2)
        qlo, qhi = q.min(0), q.max(0)
        corners = np.array([[qlo[0], qlo[1]], [qhi[0], qlo[1]], [qhi[0], qhi[1]], [qlo[0], qhi[1]]])
        frm = (o, uv, qlo, qhi, corners, float(np.linalg.norm(qhi - qlo)) + 10.0, box)

        # material on planes parallel to the face
        rot = np.array([u, v, n])
        xf = np.hstack([rot, (-rot @ o)[:, None]])
        moved = [(parts[j], parts[j].man.transform(xf)) for j in near]

        def material(height):
            """-> ({kind: owned section}, any concrete, other construction) at this height."""
            sec = {}
            for part, tm in moved:
                c = tm.slice(height)
                if not c.is_empty():
                    key = part.kind if part.cls == POURED else (-1 if part.cls == OTHER else -2)
                    sec.setdefault(key, []).append(c)
            cut_o = cs_union(sec.get(-1, []))
            cut_v = cs_union(sec.get(-2, []))
            cut = cut_o + cut_v
            kinds = {}
            prev = CrossSection()
            for k in range(nk):
                if k not in sec:
                    continue
                pk = cs_union(sec[k])
                own = pk - prev - cut
                if not own.is_empty():
                    kinds[k] = own
                prev = prev + pk
            return kinds, prev - cut, cut_o - cut_v

        # 1 mm resolution: kind PROBE behind the face, no concrete PROBE in front of it
        mat = {-1.0: material(-probe), 1.0: material(probe)}
        seen = {}
        for sign in (1.0, -1.0):
            nz = sign * n[2]
            col = 0 if nz > 0.5 else (2 if nz < -0.5 else 1)
            for k, own in mat[-sign][0].items():
                a = patch ^ own
                if a.is_empty():
                    continue
                e = a - mat[sign][1]                     # concrete in front: internal joint
                if not e.is_empty():
                    seen.setdefault(sign, []).append((k, col, e, mat[sign][2]))
        if not seen:
            continue
        # the face itself: concrete right behind it and none right in front of
        # it means it ends the concrete and stands for the boundary of its
        # stack; concrete on both sides means it is no boundary at all (a cut
        # through a part), whatever the probes meet further away
        close = {-1.0: material(-fine)[1], 1.0: material(fine)[1]}
        inside = close[-1.0] ^ close[1.0]
        ends[fi] = {}
        for sign, regions in seen.items():
            end = (patch ^ close[-sign]) - close[sign]
            if not end.is_empty():
                ends[fi][sign] = [o + poly @ uv.T for poly in end.to_polygons()]
            taken = None
            for k, col, e, other in regions:
                mine = e ^ end
                if not mine.is_empty():
                    if taken is None:
                        taken = CrossSection(stood_for(fi, sign, frm, ends, True), FillRule.Positive)
                    mine = mine - taken
                    if not mine.is_empty():
                        book(k, col, mine, other, fi, sign)
                rest = e - end - inside
                if not rest.is_empty():
                    late.append((fi, sign, frm, k, col, rest, other))

    # boundary that no face ending the concrete stands for (layers thinner than FINE)
    spare = [None] * len(faces)
    for fi, sign, frm, k, col, rest, other in late:
        o, uv = frm[0], frm[1]
        done = stood_for(fi, sign, frm, ends, False) + stood_for(fi, sign, frm, spare, True)
        left = rest - CrossSection(done, FillRule.Positive) if done else rest
        if not left.is_empty():
            book(k, col, left, other, fi, sign)
            if spare[fi] is None:
                spare[fi] = {}
            spare[fi].setdefault(sign, []).extend(o + poly @ uv.T for poly in left.to_polygons())
    return res


def takeoff(path, probe=PROBE):
    t0 = time.time()
    parts = load_parts(path)
    log('geometry: %.1f s, %d solids' % (time.time() - t0, len(parts)))
    t1 = time.time()
    vol = volumes(parts) if parts else [0.0] * len(KINDS)
    log('volumes: %.1f s' % (time.time() - t1))
    t1 = time.time()
    area = surface(parts, probe) if parts else np.zeros((len(KINDS), 4))
    log('surface: %.1f s' % (time.time() - t1))
    kinds = {}
    for k, name in enumerate(KINDS):
        kinds[name] = {
            'volume_m3': round(vol[k] / 1e9, 6) + 0.0,
            'top_m2': round(area[k, 0] / 1e6, 6) + 0.0,
            'side_m2': round(area[k, 1] / 1e6, 6) + 0.0,
            'soffit_m2': round(area[k, 2] / 1e6, 6) + 0.0,
            'contact_m2': round(area[k, 3] / 1e6, 6) + 0.0,
        }
    return {'kinds': kinds}


def main(argv):
    if len(argv) != 3:
        log('usage: pour_takeoff.py <model.ifc> <result.json>')
        return 2
    result = takeoff(argv[1])
    with open(argv[2], 'w') as fh:
        json.dump(result, fh, indent=2)
        fh.write('\n')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
POUR_TAKEOFF_PY

/opt/conda/bin/python "$OUT/pour_takeoff.py" /app/inputs/a8.ifc "$OUT/a8_takeoff.json"
