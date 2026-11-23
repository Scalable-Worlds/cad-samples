"""World-space clearance, blocking contact and motion measurements."""
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path

import FreeCAD as App
import Part


LINEAR = 0.05
ANGULAR = 0.01
CONTACT_GAP = 1e-6
ANGLE_STEP = 1e-6


def _frames(folder):
    rows = json.loads((Path(folder) / "frames.json").read_text())
    for row in rows:
        for key in ("placement", "effective"):
            coords = row[key]
            if len(coords) != 7 or not all(map(math.isfinite, coords)):
                raise ValueError("Invalid occurrence frame coordinates")
            if not any(coords[3:]):
                raise ValueError("The occurrence frame has no rotation")
            row[key] = App.Placement(App.Vector(*coords[:3]), App.Rotation(*coords[3:]))

    return rows


def _ancestry(rows, idx):
    while idx is not None:
        if not isinstance(idx, int) or not 0 <= idx < len(rows):
            raise ValueError("Invalid occurrence ancestry")
        yield idx
        parent = rows[idx]["parent"]
        if parent is not None and (not isinstance(parent, int) or not 0 <= parent < idx):
            raise ValueError("Invalid occurrence parent")
        idx = parent


def _native_roles(rows):
    # Labels identify only parts in the supplied comparator. The revised model
    # is matched by installed material and actual frames on its occurrence chain.
    def locate(*labels):
        for idx in range(len(rows)):
            chain = [rows[i]["label"] for i in _ancestry(rows, idx)]
            if chain[:len(labels)] == list(labels):
                return idx
        raise ValueError(f"Could not locate supplied native occurrence {labels}")
    return dict(
        boom=locate("Boom_Pipe", "Assembly_BoomVane"),
        vane=locate("Vane", "Assembly_BoomVane"),
        stop=locate("Tail_Stop_HighEnd", "Assembly_BoomVane"),
        pipe=locate("Pipe", "YawBearing"),
        back=locate("Blade_Assembly_BackDisk"),
        hinge_pipe=locate("Pipe", "Hinge_Outer"),
        hinge=locate("Hinge_Outer", "Tail"),
        tail=locate("Tail", "Tail_Assembly"),
    )


def _error(a, b):
    pos = (a.Base - b.Base).Length
    rot = a.Rotation.inverted() * b.Rotation
    angle = math.degrees(rot.Angle)
    return pos, min(angle, abs(360 - angle))


def _mounted_component(expected, actual, order, index, *, partial=False):
    """Pair the installed material and select its innermost mounting instance."""
    from occurrences import mounted_index
    source = expected[index]
    members = {order[i] for i in source["members"]}
    match = mounted_index(actual, members)
    if match is None and partial:
        # Independently controlled parts may be outside the tail's assembly
        # node. Their rigid motion is measured in the same native trial.
        match = mounted_index(actual, members, subset=True)
    if match is None:
        raise ValueError(f"No installed occurrence corresponds to {source['label']}")
    return match


def resolve_trial(expected_dir, actual_dir):
    """Select the loaded joint by material, and repeat a trial if its roles differ."""
    from occurrences import _pair_indices
    order = _pair_indices(expected_dir, actual_dir)
    if not order:
        raise ValueError("No one-to-one correspondence for the moving tail")
    expected, actual = _frames(expected_dir), _frames(actual_dir)
    roles = _native_roles(expected)
    tail, hinge = (_mounted_component(expected, actual, order, roles[key], partial=(key == "tail")) for key in ("tail", "hinge"))
    return None if actual[tail].get("trial_hinge") == hinge else (tail, hinge)


def _measured_members(rows, turn, solids, order, directory, required=None):
    """Compare each advanced BRep to the rigid image of its own rest material."""
    from boundary import boundary_error
    leaves = {i: row for row in rows if row["leaf"] for i in row["members"]}
    advanced = None
    def agrees(row, at, solid, moving):
        nonlocal advanced
        if len(row.get("advanced_ids", [])) != len(row["members"]):
            return False
        if required is None:
            # The comparator's own joint provides the fixed/moving class even
            # for a rotationally symmetric member centered on the hinge line.
            after_frame = App.Rotation(*row["advanced_frame"][3:])
            base = row["effective"].Rotation
            expected = turn.Rotation * base if moving else base
            angle = math.degrees((after_frame.inverted() * expected).Angle)
            if min(angle, abs(360-angle)) > ANGULAR:
                return False
        bound = row["rigid_bounds"][at][moving]
        if bound is not None and bound <= LINEAR:
            return True
        if advanced is None:
            advanced = Part.read(str(Path(directory) / "advanced.brep")).Solids
        after = advanced[row["advanced_ids"][at]]
        before = solid.copy()
        if moving:
            before.Placement = turn * before.Placement
        a, b = before.optimalBoundingBox(False), after.optimalBoundingBox(False)
        if any(abs(getattr(a, name)-getattr(b,name)) > LINEAR for name in
               ("XMin","YMin","ZMin","XMax","YMax","ZMax")):
            return False
        return boundary_error(before, after) <= LINEAR
    states = []
    for pos, (idx, solid) in enumerate(zip(order, solids, strict=True)):
        row = leaves[idx]
        at = row["members"].index(idx)
        states.append(tuple((required is None or required[pos][moving]) and agrees(row,at,solid,moving)
                            for moving in (False, True)))
    return states


def _turn(hinge):
    axis = hinge.Rotation.multVec(App.Vector(0, 0, 1))
    return App.Placement(App.Vector(), App.Rotation(axis, 20), hinge.Base)


def _world_geometry(expected, actual, pairs, order, expected_dir, actual_dir):
    roles = _native_roles(expected)
    refs, copies = pairs
    assert refs and len(refs) == len(copies), "No one-to-one installed material correspondence"
    def part(role, shapes):
        indices = expected[roles[role]]["members"]
        return shapes[indices[0]] if len(indices) == 1 else Part.makeCompound([shapes[i] for i in indices])
    indices = {role: _mounted_component(expected, actual, order, roles[role], partial=(role == "tail"))
               for role in ("back", "tail", "hinge", "stop", "vane", "pipe", "boom")}
    mounted = {role: actual[i] for role, i in indices.items()}
    parent = mounted["tail"]["placement"]
    tail = mounted["tail"]["effective"]
    hinge = mounted["hinge"]["effective"]
    back = mounted["back"]["effective"]
    frames = parent, parent.inverse() * tail, tail.inverse() * hinge
    errors = [_error(mounted[role][key], expected[roles[role]][key])
              for role, key in (("tail", "placement"), ("tail", "effective"),
                                ("hinge", "effective"), ("back", "effective"))]
    reference_hinge = expected[roles["hinge"]]["effective"]
    axis = reference_hinge.Rotation.multVec(App.Vector(0,0,1))
    actual_axis = hinge.Rotation.multVec(App.Vector(0,0,1))
    offset = hinge.Base - reference_hinge.Base
    normal = axis.cross(actual_axis)
    axis_distance = offset.cross(axis).Length
    if normal.Length:
        axis_distance = min(axis_distance, abs(offset.dot(normal)) / normal.Length)
    if "mobility" in expected[roles["tail"]]:
        reference_states = expected[roles["tail"]]["mobility"]
    else:
        reference_states = _measured_members(expected, _turn(reference_hinge), refs, range(len(refs)), expected_dir)
    measured = _measured_members(actual, _turn(hinge), copies, order, actual_dir, reference_states)
    follows = (mounted["tail"].get("trial_hinge") == indices["hinge"] and
               all(any(required and observed for required, observed in zip(r, m))
                   for r, m in zip(reference_states, measured, strict=True)))
    return (part("vane", copies), part("stop", copies), part("pipe", copies), part("pipe", refs),
            back, hinge, frames, errors, follows, axis_distance, mounted["boom"]["effective"])


def _speed_bound(shape, origin, axis):
    # Convexity gives an upper bound on every point's radius from the hinge.
    bb = shape.BoundBox
    points = product((bb.XMin, bb.XMax), (bb.YMin, bb.YMax), (bb.ZMin, bb.ZMax))
    return math.radians(max((App.Vector(*p) - origin).cross(axis).Length for p in points))


def _rotated(shape, origin, axis, angle):
    s = shape.copy()
    s.rotate(origin, axis, angle)
    return s


def _first_contact(gap, speed):
    # For a rigidly rotated part, separation is Lipschitz with constant speed.
    # Search intervals from left to right. The two endpoint cones rule out
    # intervals whose interiors cannot reach the prescribed gap.
    def search(left, right, dl, dr):
        if dl <= CONTACT_GAP:
            return left
        if (dl + dr - speed * (right - left)) / 2 > CONTACT_GAP:
            return None
        if right - left <= ANGLE_STEP:
            return right if dr <= CONTACT_GAP else None
        mid = (left + right) / 2
        dm = gap(mid)
        first = search(left, mid, dl, dm)
        return first if first is not None else search(mid, right, dm, dr)
    return search(0., 130., gap(0.), gap(130.))


def _stroke_clearance(clearance, support_gap, speed, contact, lower):
    # A separating plane supplies a cheap lower bound on Euclidean clearance.
    # Use exact BRep distances where its projection cannot certify a sample.
    # Endpoint cones bound the distance between samples over the entire stroke.
    observed = clearance(contact)
    @lru_cache(None)
    def lower_at(t):
        nonlocal observed
        bound = support_gap(t)
        if bound >= lower:
            return bound
        distance = clearance(t)
        observed = min(distance, observed)
        return distance

    def refine(left, right, dl, dr):
        if observed < lower or (dl + dr - speed * (right - left)) / 2 >= lower:
            return
        if right - left <= ANGLE_STEP:
            return
        mid = (left + right) / 2
        dm = lower_at(mid)
        refine(left, mid, dl, dm)
        refine(mid, right, dm, dr)

    refine(0., contact, lower_at(0.), observed)
    return observed


def _separation_bound(vane, origin, axis, disk_origin, disk_axis, radius, witness, travel):
    # Support of a filled circular disk is analytic. For the other side of
    # the separating plane, OCCT bounds the whole trimmed, placed vane.
    normal = (witness[0] - witness[1]).normalize()
    to_z = App.Rotation(normal, App.Vector(0, 0, 1))
    turn = App.Placement(App.Vector(), to_z)
    limit = normal.dot(disk_origin) + radius * normal.cross(disk_axis).Length

    def bound(t):
        at = _rotated(vane, origin, axis, t)
        at.Placement = turn * at.Placement
        translation = App.Rotation(axis, t).multVec(travel)
        return at.optimalBoundingBox(False, True).ZMin - limit - max(0., normal.dot(translation))
    return bound


def _wall_witness(face, centre, axis, radius, point):
    """Measure the actual trimmed contact and a two-dimensional neighbourhood."""
    surf = face.Surface
    u, v = surf.parameter(point)
    u0, u1, v0, v1 = face.ParameterRange
    def in_chart(t, lo, hi, period):
        if period:
            t += round(((lo+hi)*.5-t)/period)*period
        return min(hi,max(lo,t))
    u = in_chart(u,u0,u1,surf.UPeriod() if surf.isUPeriodic() else None)
    v = in_chart(v,v0,v1,surf.VPeriod() if surf.isVPeriodic() else None)
    if (face.valueAt(u,v)-point).Length > CONTACT_GAP:
        return False

    # A tangent plane departs by the contact-gap tolerance over a chord
    # of length sqrt(8 R gap). Sample outward from both sides of seams.
    chord = math.sqrt(8*radius*CONTACT_GAP)
    def steps(at,lo,hi,speed):
        extent = min((hi-lo)*.5, chord/speed) if speed else (hi-lo)*.5
        return {min(hi,max(lo,at+k*extent)) for k in (-1,-.5,0,.5,1)}
    us = steps(u,u0,u1,surf.getDN(u,v,1,0).Length)
    vs = steps(v,v0,v1,surf.getDN(u,v,0,1).Length)

    def radial_at(x,y):
        delta = face.valueAt(x,y)-centre
        radial = delta-axis*delta.dot(axis)
        if abs(radial.Length-radius) > CONTACT_GAP:
            return False
        radial.normalize()
        normal = face.normalAt(x,y)
        return normal.dot(radial) > 0 and radius*normal.cross(radial).Length <= CONTACT_GAP

    if not radial_at(u,v):
        return False
    for x,y in product(us,vs):
        if face.isPartOfDomain(x,y) and not radial_at(x,y):
            return False
    # A very narrow trim can truncate every off-contact transverse sample.
    # Test the local bending there as well, independent of the OCCT class.
    curvatures = [abs(k) for k in face.curvatureAt(u,v)]
    return min(curvatures) < 1/radius and max(curvatures) > 0


def _installed_wall_gap(stop, pipe, reference):
    nominal = max((f.Surface for f in reference.Faces if isinstance(f.Surface, Part.Cylinder)),
                  key=lambda s: s.Radius)
    contacts = []
    for face in pipe.Faces:
        gap, points, _ = stop.distToShape(face)
        if gap <= CONTACT_GAP:
            contacts.extend((gap,face,end) for _,end in points)

    def distances(centre,axis,radius):
        return [d for d,f,p in contacts if _wall_witness(f,centre,axis,radius,p)]
    matches = distances(nominal.Center,nominal.Axis,nominal.Radius)
    if matches:
        return min(matches)
    # Preservation permits small changes to the actual tube's radius and pose.
    # Obtain its cylinder independently from the full lateral surface or the
    # tube's circular end rims. A fitted support is only a hypothesis until
    # the actual contacting face and its neighbourhood pass the same test.
    from OCP.ShapeAnalysis import ShapeAnalysis_CanonicalRecognition
    from OCP.gp import gp_Cylinder
    from boundary import _adapt
    cylinders = [(e.Curve.Center,e.Curve.Axis,e.Curve.Radius) for e in pipe.Edges
                 if isinstance(e.Curve,Part.Circle)]
    for face in pipe.Faces:
        cyl = gp_Cylinder()
        recognizer = ShapeAnalysis_CanonicalRecognition(_adapt(face).Face())
        if recognizer.IsCylinder(CONTACT_GAP,cyl):
            cylinders.append((App.Vector(*cyl.Location().Coord()),
                              App.Vector(*cyl.Axis().Direction().Coord()), cyl.Radius()))
    for centre,axis,radius in cylinders:
        if abs(radius-nominal.Radius) > LINEAR:
            continue
        if (centre-nominal.Center).cross(nominal.Axis).Length > LINEAR:
            continue
        matches = distances(centre,axis,radius)
        if matches:
            return min(matches)
    return math.inf


def _export_errors(transforms, frames):
    # The public transform sequence has a world parent followed by the local
    # tail and hinge frames. Compare rotations as operators (unit quaternions).
    items = transforms["transforms"]
    if [item.get("name") for item in items] != ["parent", "tail", "hinge"]:
        return [(math.inf, math.inf)]
    errors = []
    for item, frame in zip(items, frames):
        rot = App.Rotation(App.Vector(*item["axis"]), math.degrees(item["angle"]))
        pose = App.Placement(App.Vector(*item["position"]), rot)
        errors.append(_error(pose, frame))
    return errors


def _disk_sweep(disk, travel, axis, radius):
    """Union of the blade disks in the coordinate system of the as-built vane."""
    if travel.Length == 0:
        return disk
    height = abs(axis.dot(travel))
    if height:
        prism = disk.extrude(travel)
        # A closed vane thicker than this entire slab has the same distance
        # to its trimmed boundary. This avoids unstable solid classification.
        return prism if height > CONTACT_GAP else Part.makeCompound(prism.Faces)
    # A translation in the disk plane sweeps a stadium, a face rather than a solid.
    in_plane = travel - axis * axis.dot(travel)
    along = in_plane / in_plane.Length
    side = axis.cross(along) * radius
    start = disk.CenterOfMass
    finish = start + in_plane
    boundary = [Part.makeLine(start+side, finish+side),
                Part.Arc(finish+side, finish+along*radius, finish-side).toShape(),
                Part.makeLine(finish-side, start-side),
                Part.Arc(start-side, start-along*radius, start+side).toShape()]
    return Part.Face(Part.Wire(boundary))


def analyze(expected_dir, actual_dir, radius, q, motion=None, pairs=None, boom_change=0.):
    if pairs is None:
        from occurrences import paired_solids
        pairs = paired_solids(expected_dir, actual_dir)
    expected = _frames(expected_dir)
    actual = _frames(actual_dir)
    from occurrences import _pair_indices
    vane, stop, pipe, ref_pipe, back, hinge, frames, errors, follows, axis_distance, boom = _world_geometry(
        expected, actual, pairs, _pair_indices(expected_dir, actual_dir), expected_dir, actual_dir)
    origin = hinge.Base
    axis = hinge.Rotation.multVec(App.Vector(0, 0, 1))
    disk = Part.Face(Part.Wire(Part.makeCircle(
        radius, back.Base, back.Rotation.multVec(App.Vector(0, 0, 1)))))
    disk_axis = back.Rotation.multVec(App.Vector(0, 0, 1))
    boom_travel = boom.Rotation.multVec(App.Vector(0, 0, boom_change))
    def stop_at(angle):
        return _rotated(stop, origin, axis, angle)
    @lru_cache(None)
    def vane_distance(angle):
        travel = App.Rotation(axis, angle).multVec(boom_travel)
        envelope = _disk_sweep(disk, travel, disk_axis, radius)
        return _rotated(vane, origin, axis, angle).distToShape(envelope)
    def clearance(angle):
        return vane_distance(angle)[0]
    gap = lambda a: stop_at(a).distToShape(pipe)[0]
    contact = _first_contact(gap, _speed_bound(stop, origin, axis))
    export_errors = _export_errors(motion, frames) if motion is not None else [(0., 0.)]
    result = dict(kinematic=follows, contact=contact, clearance=math.nan, minimum=math.nan, overlap=0.,
                  outer_gap=math.inf, wall_contact=False, datum_linear=axis_distance,
                  datum_angular=max(a for _, a in errors),
                  motion_linear=max(d for d, _ in export_errors),
                  motion_angular=max(a for _, a in export_errors))
    if contact is not None:
        at_contact = stop_at(contact)
        outer_gap = _installed_wall_gap(at_contact, pipe, ref_pipe)
        wall_contact = outer_gap <= CONTACT_GAP
        distance, witnesses, _ = vane_distance(contact)
        if distance > 0:
            bound = _separation_bound(vane, origin, axis, back.Base,
                                      disk_axis, radius, witnesses[0], boom_travel)
        else:
            bound = lambda t: 0.
        result.update(clearance=distance,
                      minimum=_stroke_clearance(clearance, bound, _speed_bound(vane, origin, axis) + math.radians(boom_travel.Length),
                                                contact, q * radius - LINEAR),
                      overlap=stop_at(contact + .1).common(pipe).Volume,
                      outer_gap=outer_gap, wall_contact=wall_contact)
    return result
