"""Read the settled native assembly and exercise its installed moving joint."""
from collections import defaultdict
from functools import cache
import json
import math
import os
from pathlib import Path
import re
import sys
from xml.etree import ElementTree
from zipfile import ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from solids import FreeCAD as App, snapshot
from frames import Inventory, mounted_index, physical_roots
from native_fields import read_native_data


def dependencies(wind):
    """FreeCAD returns the dependency closure in leaves-first recompute order."""
    return list(wind.getDependentDocuments(True))


IDENTIFIER = r"(?:<<[^>]+>>|[A-Za-z_][\w]*|\$?[A-Z]+\$?\d+)"
REFERENCE = re.compile(rf"(?:(?P<doc>{IDENTIFIER})#)?(?:(?P<obj>{IDENTIFIER})\.)?\.?(?P<prop>{IDENTIFIER})")


def _identifier(value):
    return (value or "").removeprefix("<<").removesuffix(">>").replace("$", "")


def governing_settings(wind, returned_master=None):
    """Trace model-driving expressions, including individual spreadsheet cells."""
    docs = dependencies(wind)
    documents = {key: d for d in docs for key in (d.Label, d.Name)}
    named = {d: {key: obj for obj in d.Objects for key in (obj.Label, obj.Name)} for d in docs}
    props, formulas, aliases = {}, {}, {}
    for doc in docs:
        for obj in doc.Objects:
            props[obj] = set(obj.PropertiesList)
            aliases[obj] = {}
            if obj.isDerivedFrom("Spreadsheet::Sheet"):
                for cell in obj.getNonEmptyCells():
                    name = obj.getAlias(cell) or cell
                    aliases[obj][cell] = name
                    props[obj].add(name)
                    content = obj.getContents(cell)
                    if content.startswith("="):
                        formulas[obj, name] = content[1:]
            else:
                for attr, expr in obj.ExpressionEngine:
                    formulas[obj, attr.lstrip(".")] = expr

    def references(obj, expression):
        for m in REFERENCE.finditer(expression):
            doc = documents.get(_identifier(m["doc"])) if m["doc"] else obj.Document
            other = named.get(doc, {}).get(_identifier(m["obj"])) if m["obj"] else obj
            prop = _identifier(m["prop"])
            if other is not None and prop in props.get(other, ()):
                yield other, aliases[other].get(prop, prop)
            elif m["obj"] and not m["doc"] and _identifier(m["obj"]) in props[obj]:
                # Local vector/placement property followed by a sub-property.
                base = _identifier(m["obj"])
                yield obj, aliases[obj].get(base, base)

    @cache
    def parents(node):
        return tuple(references(node[0], formulas.get(node, "")))

    def closure(seeds):
        found = set()
        pending = list(seeds)
        while pending:
            node = pending.pop()
            if node not in found:
                found.add(node)
                pending.extend(parents(node))
        return found

    seeds = set()
    installed = set()
    for root in physical_roots(wind):
        installed.update([root, *root.OutListRecursive])
    for obj, key in formulas:
        if obj not in installed or obj.isDerivedFrom("Spreadsheet::Sheet") or not hasattr(obj, "Shape"):
            continue
        prop = re.split(r"[.\[]", key, 1)[0]
        # Stock and pose depend on native modeling properties. An added shop
        # readout can itself have an expression without driving the part.
        if prop in obj.PropertiesList and 21 in obj.getPropertyStatus(prop) and not getattr(obj, "Proxy", None):
            continue
        seeds.update(parents((obj, key)))
    used = closure(seeds)

    @cache
    def scalar(node, kind):
        obj, prop = node
        raw = getattr(obj, prop, None)
        if hasattr(raw, "Unit"):
            if raw.Unit.Type not in (kind, "1"):
                return None
            raw = raw.Value
        return float(raw) if isinstance(raw, (int, float)) and math.isfinite(raw) else None

    def words(value):
        return re.sub("[^a-z0-9]", "", value.lower())

    # Start from material-driving model expressions, not every shop readout
    # that happens to carry a familiar name. In particular a multi-document
    # assembly may put exactly the same readout on many components.
    roles = {"angle": ("Angle", "maximumfurlangle"), "width": ("Length", "highendstopwidth")}
    published = {role: {(o, prop) for o in props for prop in props[o] if words(prop) == name}
                 for role, (_, name) in roles.items()}
    rotations = {}
    for node in used:
        obj, prop = node
        match = re.fullmatch(r"\s*rotation\s*\((.*)[;,]\s*([^;,]+)\)\s*",
                             formulas.get(node, ""), re.I)
        if match:
            operands = {n for n in references(obj, match[2]) if scalar(n, "Angle") is not None}
            if operands:
                rotations[node] = (operands, match[1])
    # A public unrounded readback may forward an internal angle driver. Follow
    # its identity even when neither the driver nor rotation property uses the
    # conventional spelling.
    declared_angles = closure(published["angle"])
    rotation_nodes = {n for n, (_, axis) in rotations.items()
                      if "furl" in words(n[1]) or "furl" in words(axis)}
    if not rotation_nodes:
        rotation_nodes = {n for n, (operands, _) in rotations.items() if operands & declared_angles}
    angle_inputs = {driver for n in rotation_nodes for driver in rotations[n][0]}
    width_inputs = set()
    # The supplied coped stock is a native profile. Find the installed cut
    # through its actual Base dependency. Of the high and low stop profiles,
    # the high stop's remaining cut dimensions depend on the full-furl pose.
    for cut in installed:
        if not cut.isDerivedFrom("Part::Cut") or cut.Base is None:
            continue
        for sketch in (cut.Base, *cut.Base.OutListRecursive):
            if not sketch.isDerivedFrom("Sketcher::SketchObject"):
                continue
            dims = [(key, expr) for key, expr in sketch.ExpressionEngine if key.startswith("Constraints[")]
            operands = {n for _, expr in dims for n in references(sketch, expr)}
            if not closure(operands) & rotation_nodes:
                continue
            for key, expr in dims:
                index = int(key.split("[")[1].split("]")[0])
                if sketch.Constraints[index].Type == "DistanceY":
                    width = sketch.getDatum(index).Value
                    if width > 0 and abs(width - sketch.Shape.BoundBox.YMax) <= 0.05:
                        width_inputs.update(references(sketch, expr))

    candidates = defaultdict(dict)
    for role, inputs in (("angle", angle_inputs), ("width", width_inputs)):
        kind = roles[role][0]
        tol = 0.01 if role == "angle" else 0.05
        direct = {}
        for driver in inputs:
            value = scalar(driver, kind)
            if value is None:
                continue
            for node in closure([driver]) & used:
                stored = scalar(node, kind)
                if stored is not None and abs(stored - value) <= tol:
                    # The mechanical operand supplies the measured value,
                    # including where a forwarding readout converts units.
                    direct.setdefault(node, set()).add(driver)
                    candidates[node[0].Document].setdefault(role, set()).add((1, driver))
        for node in published[role]:
            value = scalar(node, kind)
            if value is None:
                continue
            traced = closure([node])
            matched = {driver for n in traced & direct.keys() for driver in direct[n]
                       if abs(value - scalar(driver, kind)) <= tol}
            if matched:
                candidates[node[0].Document].setdefault(role, set()).update((2, driver) for driver in matched)
            elif not inputs and any(scalar(n, kind) == value for n in traced & used):
                # Other construction styles need no sketch or stock Boolean.
                # An operative public value also provides a native interface.
                candidates[node[0].Document].setdefault(role, set()).add((1, node))

    # A native master can distribute fitted data across its own dependencies.
    # Distinguish a genuine build interface (its rotor and tail inputs enter
    # the machinery) from a bearing fixture forwarding just a pose.
    input_owners = defaultdict(set)
    for o, prop in used:
        if prop in {"RotorDiameter", "VaneLength"}:
            input_owners[o.Document].add(prop)
    for doc, inputs in input_owners.items():
        if inputs != {"RotorDiameter", "VaneLength"}:
            continue
        local = candidates[doc]
        for role in roles:
            if role not in local:
                options = [option for dep in dependencies(doc) if dep != doc
                           for option in candidates.get(dep, {}).get(role, [])]
                if options:
                    local[role] = options

    providers = {d for d, values in candidates.items() if roles.keys() <= values.keys()}
    # An auxiliary document can legitimately mirror both settings. Where
    # multiple documents expose them, prefer the cross-component build input
    # provider over a fixture that merely reads the master.
    build_interfaces = {d for d in providers if input_owners[d] == {"RotorDiameter", "VaneLength"}}
    providers = build_interfaces or providers
    if returned_master is None:
        # Pick the most upstream document which actually exposes both build
        # settings. A wrapper or instrument referring to it has no such role.
        upstream = [d for d in providers if not (set(dependencies(d)) - {d}) & providers]
        if len(upstream) != 1:
            raise ValueError("Cannot distinguish the governing stop interface of the installed assembly")
        master = upstream[0]
    else:
        master = returned_master
        if master not in docs or master not in providers:
            raise ValueError("The returned master does not supply the installed stop and furl settings")

    result = []
    for role, (kind, _) in roles.items():
        options = candidates[master][role]
        top = max(priority for priority, _ in options)
        drivers = {driver for priority, driver in options if priority == top}
        values = {scalar(driver, kind) for driver in drivers}
        tolerance = 0.01 if role == "angle" else 0.05
        if max(values) > min(values) + tolerance:
            raise ValueError(f"Operative stop {role} spans {min(values)} to {max(values)}, "
                             f"exceeding {tolerance} in {kind.lower()} units")
        if len(values) == 1:
            result.append(next(iter(values)))
            continue

        # Intermediate iterates may each drive a trial rotation in the graph.
        # The released operand depends on all such predecessors. A parallel
        # calculation can instead publish their result in the master's stop interface.
        finals = {scalar(driver, kind) for driver in drivers if drivers <= closure([driver])}
        if len(finals) == 1:
            result.append(next(iter(finals)))
            continue
        field = "MaximumFurlAngle" if role == "angle" else "HighEndStopWidth"
        holder = master.getObject("HighEndStop")
        holders = [holder] if holder is not None else master.getObjectsByLabel("HighEndStop")
        readouts = {scalar((obj, field), kind) for obj in holders} - {None}
        if len(readouts) != 1:
            raise ValueError(f"No final stop {role} driver or unambiguous published master readout")
        result.append(next(iter(readouts)))
    return master, *result


def _driver_inputs(joint, expressions):
    """Find editable native inputs of the joint placement, following cell expressions."""
    docs = list(App.listDocuments().values())
    names = {n: d for d in docs for n in (d.Name, d.Label)}
    pending, seen = [(joint, expr) for _, expr in expressions], set()
    while pending:
        owner, expr = pending.pop()
        for match in REFERENCE.finditer(expr):
            doc = names.get(_identifier(match["doc"])) if match["doc"] else owner.Document
            name = _identifier(match["obj"])
            prop = _identifier(match["prop"])
            if not doc:
                continue
            obj = (doc.getObject(name) or next(iter(doc.getObjectsByLabel(name)), None)) if name else owner
            if obj is None or (obj, prop) in seen:
                continue
            seen.add((obj, prop))
            if obj.isDerivedFrom("Spreadsheet::Sheet"):
                cell = obj.getCellFromAlias(prop) or prop
                if cell not in obj.getNonEmptyCells():
                    continue
                content = obj.getContents(cell)
                if content.startswith("="):
                    pending.append((obj, content[1:]))
                    continue
            else:
                if prop not in obj.PropertiesList:
                    continue
                upstream = [e for p,e in obj.ExpressionEngine if p.strip(".").split(".")[0] == prop]
                if upstream:
                    pending.extend((obj, e) for e in upstream)
                    continue
            yield obj, prop


def _drive_joint(wind, joint, start, goal):
    import math
    def refresh():
        for doc in dependencies(wind):
            doc.recompute()
    expressions = [(p,e) for p,e in joint.ExpressionEngine if "Placement" in p]
    if not expressions:
        joint.Placement = goal
        refresh()
        return

    rotation = goal.Rotation * start.Rotation.inverted()
    def progress():
        q = (joint.Placement.Rotation * start.Rotation.inverted()).Q
        return math.degrees(2 * math.atan2(App.Vector(*q[:3]).dot(rotation.Axis), q[3]))
    def at_goal():
        angular = math.degrees((joint.Placement.Rotation.inverted() * goal.Rotation).Angle)
        return ((joint.Placement.Base - goal.Base).Length <= .05 and
                min(angular, abs(360-angular)) <= .01)

    for obj, prop in _driver_inputs(joint, expressions):
        value = getattr(obj, prop, None)
        sheet = obj.isDerivedFrom("Spreadsheet::Sheet")
        cell = (obj.getCellFromAlias(prop) or prop) if sheet else None
        saved = obj.getContents(cell) if sheet else value
        def assign(v):
            if sheet:
                obj.set(cell, str(v))
            else:
                setattr(obj, prop, v)
            refresh()
        try:
            if isinstance(value, App.Placement):
                for v in (goal * start.inverse() * value, value * start.inverse() * goal):
                    assign(v)
                    if at_goal():
                        return
            elif isinstance(value, (int,float)) or (hasattr(value, "Unit") and
                    value.Unit.Type in ("Angle", "")):
                # Measure a scalar control's gain. A builder can publish the
                # operating angle in degrees, radians, or through a gear ratio.
                x0 = value.Value if hasattr(value, "Value") else value
                suffix = " deg" if sheet and getattr(getattr(value,"Unit",None),"Type",None) == "Angle" else ""
                x1, y0 = x0 + 1, 0.
                for _ in range(12):
                    assign(str(x1) + suffix if sheet else x1)
                    if at_goal():
                        return
                    y1 = progress()
                    if y1 == y0:
                        break
                    x0, x1, y0 = x1, x1 + (20-y1)*(x1-x0)/(y1-y0), y1
        except (TypeError, ValueError, RuntimeError, ArithmeticError):
            pass
        assign(saved)
    raise ValueError("Unable to advance the installed joint through its driving input")


def motion_trial(wind, inventory, tail, hinge, folder):
    """Move the native joint with its couplings intact and capture the advanced material."""
    from itertools import product
    import Part
    joint, mount, _ = inventory.native[tail]
    hinge_frame = inventory.native[hinge][2]
    axis = hinge_frame.Rotation.multVec(App.Vector(0, 0, 1))
    turn = App.Placement(App.Vector(), App.Rotation(axis, 20), hinge_frame.Base)
    old = joint.Placement
    _drive_joint(wind, joint, old, old * mount.inverse() * turn * mount)

    advanced = Inventory(wind)
    changed = {tuple(row["path"]): row for row in advanced.frames}
    Part.makeCompound(advanced.solids).exportBrep(str(folder / "advanced.brep"))
    for row in inventory.frames:
        if not row["leaf"] or not row["members"]:
            continue
        after = changed.get(tuple(row["path"]))
        if not after or len(after["members"]) != len(row["members"]):
            continue
        row["advanced_ids"] = after["members"]
        row["advanced_frame"] = after["effective"]
        bounds = []
        for i, j in zip(row["members"], after["members"]):
            before, moved = inventory.solids[i], advanced.solids[j]
            if not before.isPartner(moved):
                bounds.append([None, None])
                continue
            delta = moved.Placement * before.Placement.inverse()
            b = before.BoundBox
            corners = [App.Vector(*p) for p in product((b.XMin,b.XMax),(b.YMin,b.YMax),(b.ZMin,b.ZMax))]
            bounds.append([max((delta.multVec(p) - op.multVec(p)).Length for p in corners)
                           for op in (App.Placement(), turn)])
        row["rigid_bounds"] = bounds
    inventory.frames[tail]["trial_hinge"] = hinge


def recompute(wind, *, force=True):
    docs = dependencies(wind)
    for doc in docs:
        if force:
            for obj in doc.Objects:
                obj.touch()
            doc.recompute(None, True, True)
        else:
            doc.recompute()
    failures = [(doc.Label, obj.Label, obj.State) for doc in docs for obj in doc.Objects
                if not obj.isValid() or "Touched" in obj.State]
    for doc in docs:
        pending = {obj.Name for obj in doc.Objects if "Proxy" in obj.PropertiesList and obj.Proxy is None}
        if pending:
            # A FeaturePython may be used without a proxy. Require restoration
            # when the native file actually declares a Python implementation.
            with ZipFile(doc.FileName) as archive:
                root = ElementTree.fromstring(archive.read("Document.xml"))
            for obj in root.findall(".//ObjectData/Object"):
                if obj.get("name") in pending and obj.find(
                        "./Properties/Property[@name='Proxy']/Python[@module]") is not None:
                    failures.append((doc.Label, obj.get("name"), "Proxy did not restore"))
    return failures


def edit_clearance(wind, package, value):
    """Change the builder's aliased input in the already-open saved master."""
    path = (package / "documents" / "Master_of_Puppets.FCStd").resolve()
    master = next((d for d in dependencies(wind) if Path(d.FileName).resolve() == path), None)
    if master is None:
        raise ValueError("The saved assembly does not link to its master input document")
    sheet = master.getObject("Spreadsheet")
    if sheet is None or not sheet.isDerivedFrom("Spreadsheet::Sheet"):
        raise ValueError("The saved master is missing its input Spreadsheet")
    cell = sheet.getCellFromAlias("VaneClearanceFraction")
    if not cell:
        raise ValueError("Missing editable VaneClearanceFraction alias in the saved master's Spreadsheet")
    sheet.set(cell, value)
    return master


def main():
    path, folder, mode, package = sys.argv[1:5]
    folder, package = Path(folder), Path(package).resolve()
    # Trusted measurement imports above precede the shipped feature modules
    # FreeCAD must restore when it opens the native production graph.
    sys.path.append(str(package.parent))
    App.ParamGet("User parameter:BaseApp/Preferences/Document").SetBool("NoPartialLoading", True)
    wind = App.openDocument(str(Path(path).resolve()))
    edited = mode in ("edited", "edited_trial")
    trial = mode in ("trial", "edited_trial")
    returned = edit_clearance(wind, package, sys.argv[5]) if edited else None
    ordinary = {}
    if edited:
        issues = recompute(wind, force=False)
        if not trial:
            from openafpm_cad_core.furl_transform import get_furl_transform
            snapshot(wind, folder / "returned")
            master, angle, width = governing_settings(wind, returned)
            ordinary = dict(returned_angle=angle, returned_width=width,
                            returned_native_data=read_native_data(master), returned_invalid_objects=issues,
                            returned_motion=get_furl_transform(wind, master))
    if mode == "saved":
        snapshot(wind, folder / "opened")
    invalid = recompute(wind)
    graph = Inventory(wind) if trial else snapshot(wind, folder)
    if trial:
        folder.mkdir(exist_ok=True, parents=True)
        motion_trial(wind, graph, *map(int, sys.argv[-2:]), folder)
    else:
        docs = dependencies(wind)
        if mode == "assembly":
            file = (package.parent / sys.argv[5]).resolve()
            returned = next((d for d in docs if Path(d.FileName).resolve() == file), None)
            if returned is None:
                raise ValueError("The returned master is absent from the reopened assembly graph")
        master, angle, width = governing_settings(wind, returned)
        external = [d.FileName for d in docs if d.FileName and not Path(d.FileName).resolve().is_relative_to(package.parent)]
        result = dict(angle=angle, width=width, external_documents=external, invalid_objects=invalid,
                      native_data=read_native_data(master), **ordinary)
        if edited:
            result["motion"] = get_furl_transform(wind, master)
        (folder / "result.json").write_text(json.dumps(result, allow_nan=False))
        # Avoid a second load for the ordinary layout. The independent matcher
        # checks both identities by installed material and repeats any other trial.
        tail = [i for i,r in enumerate(graph.frames) if r["members"] and "Tail/Tail.FCStd" in r["files"]]
        hinge = [i for i,r in enumerate(graph.frames) if r["members"] and
                 "Tail/Hinge/Outer/Tail_Hinge_Outer.FCStd" in r["files"]]
        if len(tail) == len(hinge) == 1:
            tail_index = mounted_index(graph.frames, graph.frames[tail[0]]["members"])
            hinge_index = mounted_index(graph.frames, graph.frames[hinge[0]]["members"])
            try:
                motion_trial(wind, graph, tail_index, hinge_index, folder)
            except (ValueError, RuntimeError):
                # Optimistic identification only. The paired native comparator
                # selects the actual joint for an independent trial below.
                pass
    (folder / "frames.json").write_text(json.dumps(graph.frames))
    for doc in list(App.listDocuments(True).values())[::-1]:
        App.closeDocument(doc.Name)


if __name__ == "__main__":
    main()
