"""Read material and mounted frames along native occurrence paths."""
import os
from pathlib import Path

import FreeCAD as App
import Part


def _local(obj):
    return getattr(obj, "Placement", App.Placement())


def _pose(frame):
    return [*frame.Base, *frame.Rotation.Q]


def _contents(shape):
    """Inspect every topological child, including open material beside a solid."""
    if shape.ShapeType in ("Compound", "CompSolid"):
        for part in shape.childShapes():
            yield from _contents(part)
    else:
        yield shape


def _construction(obj):
    return (obj.TypeId in {"App::Point", "App::Line", "App::Plane"} or obj.isDerivedFrom("App::OriginFeature") or obj.isDerivedFrom("App::Origin")
            or ((obj.isDerivedFrom("PartDesign::ShapeBinder") or obj.isDerivedFrom("PartDesign::SubShapeBinder"))
                and not obj.Shape.Solids)
            or obj.isDerivedFrom("PartDesign::CoordinateSystem")
            or obj.isDerivedFrom("Sketcher::SketchObject")
            or obj.isDerivedFrom("PartDesign::Plane") or obj.isDerivedFrom("PartDesign::Line")
            or obj.isDerivedFrom("PartDesign::Point"))


def _array(obj):
    link_array = obj.isDerivedFrom("App::Link") and getattr(obj, "ElementCount", 0) > 0
    draft_array = (getattr(getattr(obj, "Proxy", None), "Type", None) in
                   {"Array", "PathArray", "PointArray"} and getattr(obj, "Base", None) is not None)
    return bool((link_array or draft_array) and not getattr(obj, "Fuse", False)
                and (getattr(obj, "PlacementList", []) or getattr(obj, "ElementList", [])))


def _linked(obj):
    if not (obj.isDerivedFrom("App::Link") or obj.isDerivedFrom("App::LinkElement")):
        return None
    target = obj.LinkedObject
    if isinstance(target, tuple):
        # A subelement link exposes the selected trimmed topology via its
        # own Shape. Following the whole source would materialize extra parts.
        return None if any(target[1:]) else target[0]
    return target


def _children(obj):
    if obj.isDerivedFrom("App::LinkGroup"):
        return obj.ElementList
    if obj.isDerivedFrom("App::Part") or obj.isDerivedFrom("App::DocumentObjectGroup"):
        return obj.Group
    return None


def physical_roots(document):
    """Incoming measurement expressions do not make an assembly a subcomponent."""
    def material(obj):
        return not _construction(obj) and (_children(obj) is not None or hasattr(obj, "Shape"))

    # Only owning groups, bodies and active arrays remove their children from
    # the visible root list. A shop record can link to any of these objects.
    owned = set()
    for obj in document.Objects:
        children = _children(obj)
        if children is not None:
            owned.update(children)
        if obj.isDerivedFrom('PartDesign::Body'):
            owned.update(obj.Group)
        if _array(obj):
            owned.update(getattr(obj, 'ElementList', ()))
    return [obj for obj in document.Objects if obj not in owned and material(obj)
            and getattr(obj, "Visibility", True)]


def mounted_index(rows, members, *, subset=False):
    """Resolve the innermost installation of this material along its instance path."""
    members = frozenset(members)
    candidates = [(i, r) for i, r in enumerate(rows) if r["members"]
                  and (frozenset(r["members"]) <= members if subset else frozenset(r["members"]) == members)]
    if not candidates:
        return None
    size = max(len(r["members"]) for _, r in candidates)
    candidates = [(i, r) for i, r in candidates if len(r["members"]) == size]
    # A modeling primitive within a part may describe its shape from another
    # corner. Installations are links, bodies and parts, or standalone features.
    mounts = [(i,r) for i,r in candidates if r["mount_kind"] != "feature"] or candidates
    return max(mounts, key=lambda item: len(item[1]["path"]))[0]


class Inventory:
    def __init__(self, document):
        self.solids = []
        self.frames = []
        self.native = []
        self.datums = []
        self.invalid = []
        self.root_folder = Path(document.FileName).resolve().parent
        for root in physical_roots(document):
            if getattr(root, "Visibility", True) and not _construction(root):
                first = len(self.solids)
                self.visit(root, _local(root), None, (), (root.Name,))
                built = self.solids[first:]
                # Use FreeCAD's assembled material, including native sublinks,
                # instance scaling and any other transforms on the occurrence.
                # The path walk supplies datums; it does not replace the BRep.
                installed = Part.getShape(root).Solids
                if len(installed) != len(built):
                    raise ValueError(f"Unresolved installed material in {root.Label}: extracted {len(built)}, native {len(installed)}")
                fitted = [s.CenterOfMass for s in built]
                measured = [s.CenterOfMass for s in installed]
                # Occurrence enumeration normally follows FreeCAD's compound
                # order. Reordered aggregates need a geometric assignment.
                if any((p-q).Length > .05 for p,q in zip(fitted, measured)):
                    from scipy.optimize import linear_sum_assignment
                    costs = [[(p-q).Length for q in measured] for p in fitted]
                    _, ids = linear_sum_assignment(costs)
                    installed = [installed[i] for i in ids]
                self.solids[first:] = installed

    def visit(self, obj, mount, parent, stack, path):
        # Reusing a source at another index or via another link is an ordinary
        # instance. Detect a cycle only on the active chain of resolved edges.
        if obj in stack:
            raise ValueError(f"Cyclic assembly edge at {'.'.join(path)}")
        stack = (*stack, obj)
        target = obj
        effective = mount
        sources = []
        files = []
        last_link = None
        while not _array(target) and _linked(target) is not None:
            link = _linked(target)
            file = Path(link.Document.FileName).resolve()
            files.append(os.path.relpath(file, self.root_folder))
            sources.append(file.stem)
            last_link = effective
            if getattr(target, "LinkTransform", False):
                effective = effective * _local(link)
            target = link
            if target in stack:
                raise ValueError(f"Cyclic link at {'.'.join(path)}")
            stack = (*stack, target)
        if _construction(target):
            return
        row = dict(parent=parent, label=obj.Label, target_label=target.Label,
                   sources=sources, files=files, path=list(path),
                   mount_kind=("link" if _linked(obj) is not None else
                               "container" if _children(obj) is not None or obj.isDerivedFrom("PartDesign::Body")
                               else "feature"), placement=_pose(mount), effective=_pose(effective),
                   datum=_pose(last_link if last_link is not None and _children(target) is None
                               and not target.isDerivedFrom("PartDesign::Body") else effective),
                   members=[], leaf=False)
        index = len(self.frames)
        self.frames.append(row)
        self.native.append((obj, mount, effective))
        start = len(self.solids)
        if _array(target):
            base = getattr(target, "Base", None) or _linked(target)
            elements = getattr(target, "ElementList", [])
            positions = target.PlacementList or [_local(element) for element in elements]
            count = len(positions)
            visible = getattr(target, "VisibilityList", [])
            for i, local in enumerate(positions):
                if len(visible) == count and not visible[i]:
                    continue
                child = elements[i] if len(elements) == count else base
                if child is None:
                    raise ValueError(f"Missing source for array index {i} at {path}")
                if not elements and target.isDerivedFrom("App::Link") and target.LinkTransform:
                    local = local * _local(child)
                self.visit(child, effective * local, index, stack, (*path, str(i)))
        elif (children := _children(target)) is not None:
            for child in children:
                visible = (target.isElementVisible(child.Name) if target.isDerivedFrom("App::LinkGroup")
                           else child.Visibility)
                if visible and not _construction(child):
                    self.visit(child, effective * _local(child), index, stack, (*path, child.Name))
        elif hasattr(target, "Shape") and not target.Shape.isNull():
            # A Body exposes its finished Tip shape. A Boolean and a mirror
            # likewise expose only their finished material, never the operands.
            shape = target.Shape.copy()
            shape.Placement = effective * _local(target).inverse() * shape.Placement
            row['leaf'] = True
            for component in _contents(shape):
                if component.ShapeType != "Solid":
                    self.invalid.append([list(path), component.ShapeType])
                    continue
                self.solids.append(component)
                self.datums.append(list(effective.Rotation.Q))
        row['members'] = list(range(start, len(self.solids)))
