"""Builds reissued structural models for the rating script.

A reissue is the supplied structural model with its roof trusses replaced.
The replacement trusses are written with the constructs measured on the
supplied trusses: an IfcElementAssembly (TRUSS) placed on the roof storey,
aggregated IfcMembers with mapped Body (extruded arbitrary closed profile)
and Axis (trimmed line) representations, typed by IfcMemberType, with the
same property sets, material, style and layer.  Profile outlines are copied
from the supplied chord and web members and scaled in their own plane.
"""
import math
import random

import numpy as np

import ifcopenshell
import ifcopenshell.api

_GUID_CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_$"


def _body_solid(member):
    for rep in member.Representation.Representations:
        if rep.RepresentationIdentifier == "Body":
            return rep.Items[0].MappingSource.MappedRepresentation.Items[0]
    raise ValueError("member without Body")


def _capture_profile(curve):
    """Plain-data copy of a profile outline (polyline or composite curve)."""
    if curve.is_a("IfcPolyline"):
        return {"kind": "polyline", "points": [tuple(p.Coordinates) for p in curve.Points]}
    segments = []
    for seg in curve.Segments:
        parent = seg.ParentCurve
        if parent.is_a("IfcPolyline"):
            segments.append({"kind": "line", "same": bool(seg.SameSense),
                             "points": [tuple(p.Coordinates) for p in parent.Points]})
        else:
            circle = parent.BasisCurve
            segments.append({"kind": "arc", "same": bool(seg.SameSense),
                             "centre": tuple(circle.Position.Location.Coordinates),
                             "radius": float(circle.Radius),
                             "trim1": float(parent.Trim1[0].wrappedValue),
                             "trim2": float(parent.Trim2[0].wrappedValue),
                             "sense": bool(parent.SenseAgreement)})
    return {"kind": "composite", "segments": segments}


def _profile_extent(profile):
    xs, ys = [], []
    if profile["kind"] == "polyline":
        pts = profile["points"]
    else:
        pts = [p for seg in profile["segments"] if seg["kind"] == "line" for p in seg["points"]]
    for x, y in pts:
        xs.append(x)
        ys.append(y)
    return max(xs) - min(xs), max(ys) - min(ys)


def _dims(text):
    return ("%.1f" % text).rstrip("0").rstrip(".")


_BASES = {}


def _stripped_base(source_path):
    """The supplied model without its trusses, as STEP text, plus what the
    replacement trusses reuse from it (entity ids and profile templates)."""
    if source_path in _BASES:
        return _BASES[source_path]
    model = ifcopenshell.open(source_path)
    assemblies = [a for a in model.by_type("IfcElementAssembly") if a.PredefinedType == "TRUSS"]
    members = [m for a in assemblies for rel in a.IsDecomposedBy for m in rel.RelatedObjects]
    chord = next(m for m in members if _body_solid(m).SweptArea.OuterCurve.is_a("IfcPolyline"))
    web = next(m for m in members if _body_solid(m).SweptArea.OuterCurve.is_a("IfcCompositeCurve"))
    templates = {}
    for key, member in (("chord", chord), ("web", web)):
        solid = _body_solid(member)
        mtype = [r.RelatingType for r in member.IsDefinedBy if r.is_a("IfcRelDefinesByType")][0]
        templates[key] = {
            "profile": _capture_profile(solid.SweptArea.OuterCurve),
            "profile_name": solid.SweptArea.ProfileName,
            "family": mtype.Name.split(":")[0],
            "predefined": mtype.PredefinedType,
        }
    first = assemblies[0]
    body_rep = [r for r in chord.Representation.Representations if r.RepresentationIdentifier == "Body"][0]
    axis_rep = [r for r in chord.Representation.Representations if r.RepresentationIdentifier == "Axis"][0]
    owner = first.OwnerHistory
    keep = {
        "owner": owner.id(),
        "storey_placement": first.ObjectPlacement.PlacementRelTo.id(),
        "container": first.ContainedInStructure[0].id(),
        "origin": first.ObjectPlacement.RelativePlacement.Location.id(),
        "body_context": body_rep.ContextOfItems.id(),
        "axis_context": axis_rep.ContextOfItems.id(),
        "operator": body_rep.Items[0].MappingTarget.id(),
        "layer": [l for l in model.get_inverse(body_rep) if l.is_a("IfcPresentationLayerAssignment")][0].id(),
        "style": [s for s in model.get_inverse(_body_solid(chord)) if s.is_a("IfcStyledItem")][0].Styles[0].id(),
        "material": [r for r in model.get_inverse(chord) if r.is_a("IfcRelAssociatesMaterial")][0].id(),
    }
    assembly_name = first.Name.rsplit(":", 1)[0]
    types = {}
    for member in members:
        for rel in member.IsDefinedBy:
            if rel.is_a("IfcRelDefinesByType"):
                types[rel.RelatingType.id()] = rel.RelatingType
    for product in members + assemblies + list(types.values()):
        ifcopenshell.api.run("root.remove_product", model, product=product)
    # the removal stamps edited relationships with a new owner history; undo that
    for history in model.by_type("IfcOwnerHistory"):
        if history.id() == keep["owner"]:
            continue
        for user in list(model.get_inverse(history)):
            for index, value in enumerate(user):
                if value == history:
                    user[index] = owner
        model.remove(history)
    _BASES[source_path] = (model.to_string(), keep, templates, assembly_name)
    return _BASES[source_path]


class Reissue(object):
    """One reissued model under construction."""

    def __init__(self, source_path, seed):
        text, keep, templates, assembly_name = _stripped_base(source_path)
        self.model = ifcopenshell.file.from_string(text)
        self.rng = random.Random(seed)
        self.np_rng = np.random.default_rng(seed)
        self.tag = 700000 + 1000 * (seed % 97)
        self.templates = templates
        self.assembly_name = assembly_name
        for name, ident in keep.items():
            setattr(self, name, self.model.by_id(ident))
        self.unit_x = self.direction((1.0, 0.0, 0.0))
        self.sweep_axis = self.direction((-1.0, 0.0, 0.0))
        self.sweep_ref = self.direction((0.0, -1.0, 0.0))
        self.sweep_dir = self.direction((0.0, 0.0, -1.0))
        self.arc_ref = self.direction((1.0, 0.0))
        self._shared = {}
        self._props = {}

    # ------------------------------------------------------------- helpers
    def guid(self):
        return self.rng.choice(_GUID_CHARS[:4]) + "".join(self.rng.choice(_GUID_CHARS) for _ in range(21))

    def next_tag(self):
        self.tag += 1
        return str(self.tag)

    def new(self, kind, *args, **kwargs):
        return self.model.create_entity(kind, *args, **kwargs)

    def point(self, coords):
        return self.new("IfcCartesianPoint", tuple(float(c) for c in coords))

    def direction(self, ratios):
        return self.new("IfcDirection", tuple(float(c) for c in ratios))

    # ----------------------------------------------------------- profiles
    def _profile(self, key, scale):
        template = self.templates[key]
        prof = template["profile"]
        width, depth = _profile_extent(prof)
        if key == "chord":
            flange = abs(prof["points"][0][1] - prof["points"][1][1])
            web_t = 2.0 * min(abs(p[0]) for p in prof["points"])
            tail = template["profile_name"].split(" ", 1)[1]
            name = "%sx%sxx%sx%s %s" % (_dims(width * scale), _dims(depth * scale), _dims(flange * scale),
                                         _dims(web_t * scale), tail)
        else:
            lines = [s for s in prof["segments"] if s["kind"] == "line"]
            flange = abs(lines[0]["points"][0][1] - lines[0]["points"][1][1])
            web_t = 2.0 * min(abs(p[0]) for s in lines for p in s["points"])
            tail = template["profile_name"].split("  ", 1)[1]
            name = "%sx%sx%sx%s  %s" % (_dims(width * scale), _dims(depth * scale), _dims(flange * scale),
                                         _dims(web_t * scale), tail)
        if prof["kind"] == "polyline":
            pts = [self.point((x * scale, y * scale)) for x, y in prof["points"][:-1]]
            curve = self.new("IfcPolyline", pts + [pts[0]])
        else:
            segments = []
            for seg in prof["segments"]:
                if seg["kind"] == "line":
                    parent = self.new("IfcPolyline", [self.point((x * scale, y * scale)) for x, y in seg["points"]])
                else:
                    centre = self.point((seg["centre"][0] * scale, seg["centre"][1] * scale))
                    circle = self.new("IfcCircle", self.new("IfcAxis2Placement2D", centre, self.arc_ref),
                                      seg["radius"] * scale)
                    parent = self.new("IfcTrimmedCurve", circle,
                                      (self.new("IfcParameterValue", seg["trim1"]),),
                                      (self.new("IfcParameterValue", seg["trim2"]),),
                                      seg["sense"], "PARAMETER")
                segments.append(self.new("IfcCompositeCurveSegment", "CONTINUOUS", seg["same"], parent))
            curve = self.new("IfcCompositeCurve", segments, False)
        return self.new("IfcArbitraryClosedProfileDef", "AREA", name, curve), name

    # -------------------------------------------------------- member types
    def _member_type(self, key, scale, half, body_start, body_depth):
        """Shared representation maps and type for one member geometry."""
        ident = (key, round(scale, 9), round(half, 6), round(body_start, 6), round(body_depth, 6))
        if ident in self._shared:
            return self._shared[ident]
        template = self.templates[key]
        profile, name = self._profile(key, scale)
        position = self.new("IfcAxis2Placement3D", self.point((body_start, 0.0, 0.0)), self.sweep_axis, self.sweep_ref)
        solid = self.new("IfcExtrudedAreaSolid", profile, position, self.sweep_dir, float(body_depth))
        self.new("IfcStyledItem", solid, (self.style,), None)
        swept = self.new("IfcShapeRepresentation", self.body_context, "Body", "SweptSolid", (solid,))
        body_map = self.new("IfcRepresentationMap",
                            self.new("IfcAxis2Placement3D", self.origin, None, None), swept)
        p1, p2 = self.point((-half, 0.0, 0.0)), self.point((half, 0.0, 0.0))
        line = self.new("IfcLine", self.point((-half, 0.0, 0.0)), self.new("IfcVector", self.unit_x, 304.8))
        trimmed = self.new("IfcTrimmedCurve", line, (p1,), (p2,), True, "CARTESIAN")
        curve3d = self.new("IfcShapeRepresentation", self.axis_context, "Axis", "Curve3D", (trimmed,))
        axis_map = self.new("IfcRepresentationMap",
                            self.new("IfcAxis2Placement3D", self.origin, None, None), curve3d)
        common = self.new("IfcPropertySet", self.guid(), self.owner, "Pset_MemberCommon", None,
                          (self._prop("LoadBearing", "IfcBoolean", True),))
        mtype = self.new("IfcMemberType", self.guid(), self.owner, "%s:%s" % (template["family"], name), None, None,
                         (common,), (body_map, axis_map), self.next_tag(), None, template["predefined"])
        rel = self.new("IfcRelDefinesByType", self.guid(), self.owner, None, None, (), mtype)
        self.layer.AssignedItems = tuple(self.layer.AssignedItems) + (swept, curve3d)
        record = {"body_map": body_map, "axis_map": axis_map, "rel": rel, "name": name,
                  "family": template["family"]}
        self._shared[ident] = record
        return record

    def _prop(self, name, kind, value):
        key = (name, kind, value)
        if key not in self._props:
            self._props[key] = self.new("IfcPropertySingleValue", name, None, self.new(kind, value), None)
        return self._props[key]

    # -------------------------------------------------------------- members
    def _member(self, parent_placement, key, scale, mid, x_dir, z_dir, half, cut_start, cut_end, slope):
        length = 2.0 * half
        if key == "chord":
            body_depth = math.floor(length * 10.0) / 10.0
            body_start = -round(body_depth / 2.0, 2)
            span = body_depth
        else:
            body_start = -half + cut_start
            body_depth = round(length - cut_start - cut_end, 1)
            span = math.floor(length * 10.0) / 10.0
        shared = self._member_type(key, scale, half, body_start, body_depth)
        placement = self.new("IfcLocalPlacement", parent_placement,
                             self.new("IfcAxis2Placement3D", self.point(mid), self.direction(z_dir),
                                      self.direction(x_dir)))
        body = self.new("IfcShapeRepresentation", self.body_context, "Body", "MappedRepresentation",
                        (self.new("IfcMappedItem", shared["body_map"], self.operator),))
        axis = self.new("IfcShapeRepresentation", self.axis_context, "Axis", "MappedRepresentation",
                        (self.new("IfcMappedItem", shared["axis_map"], self.operator),))
        self.layer.AssignedItems = tuple(self.layer.AssignedItems) + (body, axis)
        shape = self.new("IfcProductDefinitionShape", None, None, (body, axis))
        tag = self.next_tag()
        member = self.new("IfcMember", self.guid(), self.owner, "%s:%s:%s" % (shared["family"], shared["name"], tag),
                          None, None, placement, shape, tag)
        reference = self._prop("Reference", "IfcIdentifier", shared["name"])
        common = self.new("IfcPropertySet", self.guid(), self.owner, "Pset_MemberCommon", None, (
            self._prop("IsExternal", "IfcBoolean", False), self._prop("LoadBearing", "IfcBoolean", True),
            self._prop("Slope", "IfcPlaneAngleMeasure", float(slope)), reference,
            self._prop("Span", "IfcPositiveLengthMeasure", float(span))))
        takeoff = self.new("IfcPropertySet", self.guid(), self.owner, "Pset_QuantityTakeOff", None, (reference,))
        for pset in (common, takeoff):
            self.new("IfcRelDefinesByProperties", self.guid(), self.owner, None, None, (member,), pset)
        shared["rel"].RelatedObjects = tuple(shared["rel"].RelatedObjects) + (member,)
        self.material.RelatedObjects = tuple(self.material.RelatedObjects) + (member,)
        return member

    # --------------------------------------------------------------- truss
    def add_truss(self, spec):
        """Add one truss.  `spec` keys:

        start      plan position [x, y] of the local start of the chords (mm)
        azimuth    plan angle of the chord direction in degrees, or
        direction  plan direction [dx, dy] of the chords (takes precedence)
        level      height of the lower chord axis above the roof storey (mm)
        depth      distance between the chord axes (mm)
        panels     panel runs from the local start (mm)
        chords     scale factors [lower, upper]
        posts      scale factor(s) of the posts, cycled along the truss
        diagonals  scale factor(s) of the diagonals, cycled along the truss
        skew       largest out-of-plane offset of a web end from its joint (mm)
        slip       largest in-plane offset of a web end from its joint (mm)
        shift      [lower, upper]: distance along the chords, in the direction of the panel runs, from the
                   post end at a joint of that chord to every diagonal end there (mm); at an end of the
                   chords the diagonal end stands on the chord side of the post end
        order      "listed" or "mixed": order of the members in the aggregate
        """
        panels = [float(p) for p in spec["panels"]]
        depth = float(spec["depth"])
        level = float(spec.get("level", 100.0))
        if "direction" in spec:
            d = np.array([spec["direction"][0], spec["direction"][1], 0.0], dtype=float)
        else:
            az = math.radians(spec["azimuth"])
            d = np.array([math.cos(az), math.sin(az), 0.0])
        d /= np.linalg.norm(d)
        up = np.array([0.0, 0.0, 1.0])
        normal = np.cross(d, up)
        start = np.array([spec["start"][0], spec["start"][1], level], dtype=float)
        stations = np.concatenate([[0.0], np.cumsum(panels)])
        span = float(stations[-1])
        n = len(panels)
        low_scale, top_scale = [float(v) for v in spec["chords"]]

        def cycle(value, count):
            values = value if isinstance(value, (list, tuple)) else [value]
            return [float(values[i % len(values)]) for i in range(count)]

        post_scales = cycle(spec.get("posts", 1.0), n + 1)
        diag_scales = cycle(spec.get("diagonals", 1.0), n)
        skew = float(spec.get("skew", 0.0))
        slip = float(spec.get("slip", 0.0))
        shift_low, shift_top = [float(v) for v in spec.get("shift", (0.0, 0.0))]

        def shifted(value, station):
            """Shift of the diagonal ends at a station; at an end of the chords it points into them."""
            if station == 0:
                return abs(value)
            if station == n:
                return -abs(value)
            return value

        placement = self.new("IfcLocalPlacement", self.storey_placement,
                             self.new("IfcAxis2Placement3D", self.origin, None, None))
        tag = self.next_tag()
        assembly = self.new("IfcElementAssembly", self.guid(), self.owner, "%s:%s" % (self.assembly_name, tag),
                            None, None, placement, None, tag, "NOTDEFINED", "TRUSS")
        takeoff = self.new("IfcPropertySet", self.guid(), self.owner, "Pset_QuantityTakeOff", None,
                           (self._prop("Reference", "IfcIdentifier", self.assembly_name.split(":")[-1]),))
        self.new("IfcRelDefinesByProperties", self.guid(), self.owner, None, None, (assembly,), takeoff)
        self.container.RelatedElements = tuple(self.container.RelatedElements) + (assembly,)

        chord_depth = _profile_extent(self.templates["chord"]["profile"])[1]
        web_depth = _profile_extent(self.templates["web"]["profile"])[1]
        chords, posts, diagonals = [], [], []
        for scale, height in ((low_scale, 0.0), (top_scale, depth)):
            mid = start + d * span / 2.0 + up * height
            chords.append(self._member(placement, "chord", scale, mid, d, up, span / 2.0, 0.0, 0.0, 0.0))

        def wobble():
            if skew == 0.0 and slip == 0.0:
                return np.zeros(3), 0.0
            off = normal * self.np_rng.uniform(-skew, skew)
            off = off + d * self.np_rng.uniform(-slip, slip) + up * self.np_rng.uniform(-slip, slip)
            return off, float(self.np_rng.uniform(-slip, slip))

        clear_low = 0.5 * chord_depth * low_scale + 12.7
        clear_top = 0.5 * chord_depth * top_scale + 12.7
        for i in range(n + 1):
            off, stretch = wobble()
            mid = start + d * stations[i] + up * depth / 2.0 + off
            cut_start = 0.0 if i == 0 else round(clear_low, 1)
            posts.append(self._member(placement, "web", post_scales[i], mid, up, d, depth / 2.0 + stretch,
                                      cut_start, round(clear_top, 1), 0.0))
        for i in range(n):
            rising = i < n / 2.0
            if rising:
                lead, trail = shifted(shift_low, i), shifted(shift_top, i + 1)
            else:
                lead, trail = shifted(shift_top, i), shifted(shift_low, i + 1)
            run = panels[i] + (trail - lead)
            length = math.hypot(run, depth)
            cos_a, sin_a = run / length, depth / length
            if rising:
                x_dir = d * cos_a + up * sin_a
                z_dir = d * sin_a - up * cos_a
                first_clear, last_clear = clear_low, clear_top
            else:
                x_dir = d * cos_a - up * sin_a
                z_dir = -d * sin_a - up * cos_a
                first_clear, last_clear = clear_top, clear_low
            off, stretch = wobble()
            mid = start + d * (stations[i] + lead + run / 2.0) + up * depth / 2.0 + off
            lap = 0.5 * web_depth * diag_scales[i] * cos_a / sin_a
            cuts = [round(c / sin_a + lap, 1) for c in (first_clear, last_clear)]
            diagonals.append(self._member(placement, "web", diag_scales[i], mid, x_dir, z_dir, length / 2.0 + stretch,
                                          cuts[0], cuts[1], math.degrees(math.atan2(depth, run))))
        ordered = chords + posts + diagonals
        if spec.get("order", "listed") == "mixed":
            self.rng.shuffle(ordered)
        self.new("IfcRelAggregates", self.guid(), self.owner, None, None, assembly, tuple(ordered))
        return assembly

    def write(self, path):
        self.model.write(path)


def build_reissue(source_path, case, dest_path):
    """Write the reissue described by `case` (a dict with seed and trusses)."""
    reissue = Reissue(source_path, int(case["seed"]))
    for spec in case["trusses"]:
        reissue.add_truss(spec)
    reissue.write(dest_path)
    return dest_path
