#!/bin/bash
set -euo pipefail

cd /app
mkdir -p /app/out
rm -rf -- /app/out/openafpm_cad_core
cp -a /app/inputs/openafpm_cad_core /app/out/openafpm_cad_core
cd /app/out/openafpm_cad_core
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1

/usr/bin/git apply - <<'PROJECT_EDITS'
--- a/dimension_tables.py
+++ b/dimension_tables.py
@@ -6,8 +6,7 @@
 
 from .find_descendent_by_label import find_descendent_by_label
 from .find_object_by_label import find_object_by_label
-from .load import load_alernator
-from .load_spreadsheet_document import load_spreadsheet_document
+from .load import load_turbine
 from .parameter_groups import FurlingParameters, MagnafpmParameters, UserParameters
 from .wind_turbine_shape import WindTurbineShape
 
@@ -36,10 +35,11 @@
     user_parameters: UserParameters,
     img_path_prefix: str = "",
 ) -> List[Element]:
-    spreadsheet_document = load_spreadsheet_document(
+    turbine_document, spreadsheet_document = load_turbine(
         magnafpm_parameters, furling_parameters, user_parameters
     )
-    alternator_document = load_alernator(recompute_all=True)
+    alternator_link = find_object_by_label(turbine_document, "Alternator")
+    alternator_document = alternator_link.LinkedObject.Document
     return get_dimension_tables(
         spreadsheet_document, alternator_document, img_path_prefix
     )
--- a/furl_transform.py
+++ b/furl_transform.py
@@ -1,10 +1,8 @@
-from pathlib import Path
 from typing import List, TypedDict
 
-import FreeCAD as App
-from FreeCAD import Console, Document, Placement
+from FreeCAD import Document, Placement
 
-from .find_object_by_label import find_object_by_label
+from .furl_geometry import HINGE_PATH, TAIL_PATH
 from .load import load_turbine
 from .parameter_groups import FurlingParameters, MagnafpmParameters, UserParameters
 
@@ -66,23 +64,16 @@
 
 
 def get_furl_transforms(root_document: Document) -> List[Transform]:
-    root_document_path = Path(root_document.FileName)
-    documents_path = root_document_path.parent
-    tail_document_path = documents_path.joinpath("Tail", "Tail.FCStd")
-    tail_document = App.openDocument(str(tail_document_path))
-    tail = find_object_by_label(tail_document, "Tail")
-    if len(tail.InList) == 0:
-        Console.PrintWarning(f"{tail.Label} has no parents.\n")
-        return None
-    if len(tail.InList) > 1:
-        Console.PrintWarning(f"{tail.Label} has more than 1 parent. Choosing 1st.\n")
-    tail_parent = tail.InList[0]
-    parent_placement = calculate_global_placement(tail_parent)
-    hinge_outer = find_object_by_label(tail_document, "Hinge_Outer")
+    assembly = root_document.Part
+    # Keep the published parent/tail/hinge factorization. Type 3 ends at
+    # the link placement; type 5 includes the linked object's own placement.
+    parent = assembly.getSubObject(TAIL_PATH, 3)
+    tail = assembly.getSubObject(TAIL_PATH, 5)
+    hinge = assembly.getSubObject(HINGE_PATH, 5)
     return [
-        placement_to_dict("parent", parent_placement),
-        placement_to_dict("tail", tail.Placement),
-        placement_to_dict("hinge", hinge_outer.Placement),
+        placement_to_dict("parent", parent),
+        placement_to_dict("tail", parent.inverse() * tail),
+        placement_to_dict("hinge", tail.inverse() * hinge),
     ]
 
 
@@ -93,19 +84,3 @@
         "axis": list(placement.Rotation.Axis),
         "angle": placement.Rotation.Angle,
     }
-
-
-def calculate_global_placement(child: object, placements: Placement = []) -> Placement:
-    placements.append(child.Placement)
-    in_list = child.InList
-    num_in = len(in_list)
-    if len(in_list) == 0:
-        global_placement = Placement()
-        placements.reverse()  # Reverse list in order of parent to child.
-        for placement in placements:
-            global_placement *= placement
-        return global_placement
-    if num_in > 1:
-        Console.PrintWarning(f"{child.Label} has more than 1 parent. Choosing 1st.\n")
-    parent = in_list[0]
-    return calculate_global_placement(parent, placements)
--- a/get_default_parameters.py
+++ b/get_default_parameters.py
@@ -1,4 +1,5 @@
 """Module for retrieving values for various wind turbine preset designs."""
+from copy import deepcopy
 from typing import Dict, List, Union
 
 from .wind_turbine_shape import WindTurbineShape
@@ -12,7 +13,7 @@
     key = preset if isinstance(preset, str) else preset.value
     parameters = preset_by_name[key]
     if 'inheritsFrom' not in parameters:
-        return parameters
+        return deepcopy(parameters)
     else:
         parent_parameters = preset_by_name[parameters['inheritsFrom']]
         merged_parameters = {'description': parameters['description']}
@@ -21,7 +22,7 @@
                 merged_parameters[group] = parent_parameters[group] | parameters[group]
             else:
                 merged_parameters[group] = parent_parameters[group]
-        return merged_parameters
+        return deepcopy(merged_parameters)
 
 
 def get_presets() -> List[str]:
--- a/load_root_document.py
+++ b/load_root_document.py
@@ -92,6 +92,22 @@
     if scaled_callback:
         scaled_callback("Recomputing documents", 70)
     recompute_all_documents(scaled_callback, cancel_event)
+    for document in root_documents:
+        if Path(document.FileName).resolve() == (documents_path / "WindTurbine.FCStd").resolve():
+            from .clearance_fit import fit_high_end_stop
+            from .native_stop_fit import connect_fitted_stop
+            fit_high_end_stop(
+                document, spreadsheet_document,
+                spreadsheet_document.Spreadsheet.VaneClearanceFraction,
+            )
+            groups = (
+                magnafpm_parameters,
+                dict(furling_parameters, VaneClearanceFraction=spreadsheet_document.Spreadsheet.VaneClearanceFraction),
+                user_parameters,
+            )
+            connect_fitted_stop(spreadsheet_document, groups)
+            for dependency in document.getDependentDocuments(True):
+                dependency.recompute()
 
     if cancel_event is not None and cancel_event.is_set():
         close_all_documents()
--- a/load_spreadsheet_document.py
+++ b/load_spreadsheet_document.py
@@ -1,11 +1,9 @@
-"""Module for loading spreadsheet document."""
+"""Load the master spreadsheet with the production assembly's fitted stop."""
 
 from FreeCAD import Document
 
-from .get_documents_path import get_documents_path
-from .parameter_groups import (FurlingParameters, MagnafpmParameters,
-                               UserParameters)
-from .upsert_spreadsheet_document import upsert_spreadsheet_document
+from .load import load_turbine
+from .parameter_groups import FurlingParameters, MagnafpmParameters, UserParameters
 
 __all__ = ['load_spreadsheet_document']
 
@@ -13,10 +11,5 @@
 def load_spreadsheet_document(magnafpm_parameters: MagnafpmParameters,
                               furling_parameters: FurlingParameters,
                               user_parameters: UserParameters) -> Document:
-    name = 'Master_of_Puppets'
-    documents_path = get_documents_path()
-    spreadsheet_document_path = documents_path.joinpath(f'{name}.FCStd')
-    return upsert_spreadsheet_document(spreadsheet_document_path,
-                                       magnafpm_parameters,
-                                       furling_parameters,
-                                       user_parameters)
+    _, master = load_turbine(magnafpm_parameters, furling_parameters, user_parameters)
+    return master
--- a/parameter_groups.py
+++ b/parameter_groups.py
@@ -3,7 +3,7 @@
 The wind turbine has many parameters which are organized into 3 broad categories.
 """
 
-from typing import TypedDict
+from typing import NotRequired, TypedDict
 
 __all__ = [
     "MagnafpmParameters",
@@ -209,6 +209,11 @@
     — page 60, Glossary section of "A Wind Turbine Recipe Book (2014)".
     """
 
+    VaneClearanceFraction: NotRequired[float]
+    """Least full-furl vane clearance across the interchangeable boom lengths,
+    expressed as a fraction of rotor radius. Range 0.10–0.15, default 0.10.
+    """
+
     VerticalPlaneAngle: int  # TailHingeAngle?
     """Angle between outer pipe of yaw-bearing and inner pipe of tail hinge (in degrees).
 
@@ -222,7 +227,8 @@
     """
 
     BoomLength: int  # BoomPipeLength
-    """Length of tail boom pipe.
+    """Cut length in mm. For T, H and Star builds the service family spans
+    this cut and the preset length. Cuts may be 60–140% of the preset length.
 
     See "Tail boom" section on page 31 of "A Wind Turbine Recipe Book (2014)".
     """
--- a/upsert_spreadsheet_document.py
+++ b/upsert_spreadsheet_document.py
@@ -1,3 +1,5 @@
+from decimal import Decimal
+from numbers import Real
 from pathlib import Path
 from typing import Dict, List
 
@@ -45,6 +47,10 @@
     furling_parameters: FurlingParameters,
     user_parameters: UserParameters,
 ) -> Dict[str, List[List[Cell]]]:
+    fraction = furling_parameters.get("VaneClearanceFraction", 0.10)
+    if not isinstance(fraction, (Real, Decimal)) or not 0.10 <= float(fraction) <= 0.15:
+        raise ValueError("VaneClearanceFraction must be a real number from 0.10 through 0.15")
+    furling_parameters = dict(furling_parameters, VaneClearanceFraction=float(fraction))
     rotor_disk_radius = magnafpm_parameters['RotorDiskRadius']
     calculated_wind_turbine_shape = map_rotor_disk_radius_to_wind_turbine_shape(rotor_disk_radius).to_string()
     cells = parameters_by_key_to_cells(
@@ -88,6 +94,10 @@
         raise InterruptedError("Operation was cancelled")
         
     populate_spreadsheets(document, cells_by_spreadsheet_name, cancel_event)
+    # Repopulating the master restores the stock design equations. Each new
+    # production build fits in its own context and then reconnects the driver.
+    if document.getObject("StopFit") is not None:
+        document.removeObject("StopFit")
     
     if cancel_event is not None and cancel_event.is_set():
         close_all_documents()
PROJECT_EDITS

cat > clearance_fit.py <<'PYTHON_MODULE'
"""Size the high-end stop in the assembled, resting wind turbine."""

import math

import FreeCAD as App
import Part
from scipy.optimize import root

from .furl_geometry import FurlGeometry
from .get_default_parameters import get_default_parameters
from .wind_turbine_shape import WindTurbineShape


def _rotate(shape, geometry, angle):
    moved = shape.copy()
    moved.rotate(geometry.hinge_origin, geometry.hinge_axis, angle)
    return moved


def _service_clearance(geometry, angle, length_difference):
    """Minimum over every service length between the preset and the cut.

    A translated vane against a fixed disk has the same clearance as a fixed
    vane against the oppositely translated disk. After rotating the vane and
    the boom axis, sweep the filled disk once along this relative path. Its
    end faces test the interval endpoints and its lateral face tests the
    continuum between them.
    """
    vane = _rotate(geometry.vane, geometry, angle)
    disk = geometry.disk
    if length_difference == 0:
        return vane.distToShape(disk)[0]

    rotation = App.Rotation(geometry.hinge_axis, angle)
    sweep = rotation.multVec(geometry.boom_axis) * length_difference
    disk_normal = disk.Placement.Rotation.multVec(App.Vector(0, 0, 1))
    if abs(disk_normal.dot(sweep)) <= 1e-7:
        # The zero-height limit of the oblique solid is a planar capsule.
        # Evaluate its two ends and the ruled rectangle on the two diameters.
        # In the kernel's coplanarity tolerance its distance differs from
        # the exact swept set by at most the normal displacement (1e-7 mm).
        center = disk.Placement.Base
        offset = disk_normal.cross(sweep)
        offset = offset * (disk.Edges[0].Curve.Radius / offset.Length)
        a, b = center + offset, center - offset
        rectangle = Part.Face(Part.makePolygon([a, b, b + sweep, a + sweep, a]))
        far_disk = disk.copy()
        far_disk.translate(sweep)
        service_disk = Part.makeCompound([disk, rectangle, far_disk])
    else:
        service_disk = disk.extrude(sweep)
    return vane.distToShape(service_disk)[0]


def _outer_wall(pipe):
    return max((face for face in pipe.Faces if isinstance(face.Surface, Part.Cylinder)),
               key=lambda face: face.Surface.Radius)


def _wall_residual(stop, pipe):
    # At a contact on the lateral wall, the closest radial distance equals
    # the wall radius. Unlike a solid-solid gap, this residual crosses zero.
    wall = _outer_wall(pipe)
    cylinder = wall.Surface
    _, _, z0, z1 = wall.ParameterRange
    axis = Part.makeLine(cylinder.Center + z0 * cylinder.Axis,
                         cylinder.Center + z1 * cylinder.Axis)
    return stop.distToShape(axis)[0] - cylinder.Radius


def _travel_per_degree(shape, geometry):
    # The perpendicular distance of each point from the hinge axis is bounded
    # by that of the box centre plus the box's circumscribed sphere radius.
    bounds = shape.BoundBox
    offset = bounds.Center - geometry.hinge_origin
    max_radius = offset.cross(geometry.hinge_axis).Length + 0.5 * bounds.DiagonalLength
    return math.radians(max_radius)


def _first_contact(geometry):
    """Advance by certified open clearance until first contact is resolved."""
    contact_gap = 1e-6
    angular_resolution = 1e-7
    speed = _travel_per_degree(geometry.stop, geometry)
    angle = 0.0
    while True:
        moved = _rotate(geometry.stop, geometry, angle)
        separation = moved.distToShape(geometry.pipe)[0]
        if separation <= contact_gap:
            if angle == 0:
                raise ValueError("High-end stop is in contact at rest")
            if moved.distToShape(_outer_wall(geometry.pipe))[0] > contact_gap:
                raise ValueError("First contact misses the outer wall of the yaw pipe")
            return angle
        if angle >= 130.0:
            raise ValueError("High-end stop misses the yaw pipe through 130 degrees")
        # No point of the stop can close the measured gap faster than this
        # speed. At the final approach resolve the threshold to 0.0000001 deg.
        safe_step = (separation - contact_gap) / speed
        angle = min(130.0, angle + max(safe_step, angular_resolution))


def _check_vane_stroke(geometry, contact, end_gap, target, length_difference):
    """Bound the vane clearance throughout the continuous approach to contact."""
    speed = _travel_per_degree(geometry.vane, geometry) + math.radians(abs(length_difference))
    minimum = target - 0.05

    def gap(angle):
        return _service_clearance(geometry, angle, length_difference)

    intervals = [(0.0, contact, gap(0.0), end_gap)]
    while intervals:
        lo, hi, left, right = intervals.pop()
        if min(left, right) < minimum:
            raise ValueError("Insufficient vane clearance before first stop contact")
        # Distance to a fixed disk is Lipschitz with this same speed. Use
        # both endpoint bounds and subdivide only where clearance is close.
        if (left + right - speed * (hi - lo)) / 2 >= minimum:
            continue
        mid = (lo + hi) / 2
        if mid == lo or mid == hi:
            raise ValueError("Could not resolve the vane clearance over the stroke")
        middle = gap(mid)
        intervals.extend([(lo, mid, left, middle), (mid, hi, middle, right)])


def fit_high_end_stop(wind_turbine_doc, spreadsheet_doc, clearance_fraction):
    """Set the two high-stop dimensions and propagate them to the linked parts."""
    sheet = spreadsheet_doc.HighEndStop
    inputs = spreadsheet_doc.Spreadsheet
    frame = WindTurbineShape.from_string(inputs.CalculatedWindTurbineShape)
    nominal_length = get_default_parameters(frame)["furling"]["BoomLength"]
    length_difference = float(inputs.BoomLength) - nominal_length
    radius = float(inputs.RotorDiameter) / 2
    target = clearance_fraction * radius
    angle_cell = sheet.getCellFromAlias("MaximumFurlAngle")
    width_cell = sheet.getCellFromAlias("HighEndStopWidth")
    documents = list(App.listDocuments(True).values())

    def set_dimensions(angle, width):
        sheet.set(angle_cell, f"{angle:.17g} deg")
        sheet.set(width_cell, f"{width:.17g}")
        # Yaw pipe length and assembly datums can change with the stop. Update
        # all touched documents, leaf parts before assemblies, on every trial.
        for document in documents:
            document.recompute()

    def residual(dimensions):
        angle, width = dimensions
        set_dimensions(angle, width)
        g = FurlGeometry.from_assembly(wind_turbine_doc, radius)
        clearance = _service_clearance(g, angle, length_difference)
        stop = _rotate(g.stop, g, angle)
        return clearance - target, _wall_residual(stop, g.pipe)

    initial = (sheet.MaximumFurlAngle.Value, float(sheet.HighEndStopWidth))
    solution = root(residual, initial, options={"eps": 1e-9, "xtol": 1e-10})
    angle, width = solution.x
    if (not math.isfinite(width) or width <= 0 or not 0 < angle < 130
            or max(map(abs, residual(solution.x))) > 1e-4):
        raise ValueError("Could not fit the high-end stop to the requested clearance")

    g = FurlGeometry.from_assembly(wind_turbine_doc, radius)
    contact = _first_contact(g)
    clearance = _service_clearance(g, contact, length_difference)
    overlap = _rotate(g.stop, g, contact + 0.1).common(g.pipe).Volume
    if abs(contact - angle) > 0.01 or abs(clearance - target) > 0.05 or overlap <= 0:
        raise ValueError("Fitted high-end stop failed the blocking-contact check")
    _check_vane_stroke(g, contact, clearance, target, length_difference)
PYTHON_MODULE

cat > furl_geometry.py <<'PYTHON_MODULE'
"""Placed geometry and datums for the production furling mechanism."""

from dataclasses import dataclass

import FreeCAD as App
import Part


# Use the assembly path, including the transform of each linked App::Part.
TAIL_PATH = "$Tail_Assembly.$Tail."
HINGE_PATH = TAIL_PATH + "$Hinge_Outer."
BOOM_PATH = TAIL_PATH + "$Assembly_BoomVane."


def stop_shape(wind_turbine):
    """Return the finished plate, after the clearance cut for the outer hinge."""
    return wind_turbine.Part.getSubObject(BOOM_PATH + "$Tail_Stop_HighEnd.")


@dataclass
class FurlGeometry:
    vane: Part.Shape
    stop: Part.Shape
    pipe: Part.Shape
    disk: Part.Shape
    boom_axis: App.Vector
    hinge_origin: App.Vector
    hinge_axis: App.Vector

    @classmethod
    def from_assembly(cls, wind_turbine, radius):
        assembly = wind_turbine.Part
        # Return type 5 includes the source placement for LinkTransform links.
        hinge = assembly.getSubObject(HINGE_PATH, 5)
        back_disk = assembly.getSubObject("$Alternator.$Blade_Assembly_BackDisk.", 5)
        boom_pipe = assembly.getSubObject(BOOM_PATH + "$Boom_Pipe.", 5)
        disk = Part.Face(Part.Wire(Part.makeCircle(radius)))
        disk.Placement = back_disk
        return cls(
            assembly.getSubObject(BOOM_PATH + "$Vane."),
            stop_shape(wind_turbine),
            assembly.getSubObject("$YawBearing.$Pipe."),
            disk,
            boom_pipe.Rotation.multVec(App.Vector(0, 0, 1)),
            hinge.Base,
            hinge.Rotation.multVec(App.Vector(0, 0, 1)),
        )
PYTHON_MODULE

cat > native_stop_fit.py <<'PYTHON_MODULE'
"""Keep the fitted stop parametric when a builder reopens the native project."""

import json
from pathlib import Path
import subprocess
import sys
import sysconfig

from FreeCAD import Units


def _read_inputs(sheet, groups):
    """Read the same three groups the project loader uses, from the open master."""
    def read(name):
        value = getattr(sheet, name)
        return value.Value if isinstance(value, Units.Quantity) else value

    return [{name: read(name) for name in group} for group in groups]


class StopFit:
    """An upstream FreeCAD feature. Its only dependency is the input sheet."""

    def execute(self, obj):
        previous = json.loads(obj.FittedInputs)
        current = _read_inputs(obj.Inputs, previous)
        if current == previous:
            return

        # Trial dimensions change the downstream pipe length and hinge pose.
        # Evaluate the linked model in a second kernel so a trial never writes
        # upstream into the document that FreeCAD is currently recomputing.
        python = Path(sysconfig.get_path("scripts")) / "python"
        package_parent = Path(__file__).resolve().parent.parent
        request = json.dumps(current, allow_nan=False)
        process = subprocess.run(
            [str(python), "-B", "-m", "openafpm_cad_core.native_stop_fit"],
            input=request, text=True, capture_output=True, cwd=package_parent,
        )
        if process.returncode:
            raise RuntimeError(f"Furl-stop refit failed:\n{process.stderr.strip()}")
        # FreeCAD can also print document messages to standard output.
        result = next((line for line in process.stdout.splitlines()
                       if line.startswith("OpenAFPM stop fit: ")), None)
        if result is None:
            raise RuntimeError("FreeCAD returned no fitted furl-stop dimensions")
        angle, width = json.loads(result.removeprefix("OpenAFPM stop fit: "))
        obj.Travel = angle
        obj.StockWidth = width
        obj.FittedInputs = request

    def dumps(self):
        return None

    def loads(self, state):
        pass


def connect_fitted_stop(master, groups):
    """Publish an already measured fit and retain its dependency on the inputs."""
    inputs = master.Spreadsheet
    stop = master.HighEndStop
    feature = master.addObject("App::FeaturePython", "StopFit")
    feature.Label = "Fitted high-end stop"
    feature.addProperty("App::PropertyLink", "Inputs", "Fit", "Authoritative build inputs")
    feature.addProperty("App::PropertyAngle", "Travel", "Fit", "Travel to the first blocking contact")
    feature.addProperty("App::PropertyLength", "StockWidth", "Fit", "Transverse width of the stop stock")
    # Keep the last solved input state with its two outputs in the FCStd. A
    # forced recompute of an unchanged build can reuse its verified dimensions.
    feature.addProperty("App::PropertyString", "FittedInputs", "Fit")
    feature.Inputs = inputs
    feature.Travel = stop.MaximumFurlAngle
    feature.StockWidth = stop.HighEndStopWidth
    feature.FittedInputs = json.dumps(_read_inputs(inputs, groups), allow_nan=False)
    for name in ("Inputs", "Travel", "StockWidth"):
        feature.setEditorMode(name, 1)
    feature.setEditorMode("FittedInputs", 2)
    feature.Proxy = StopFit()
    stop.set(stop.getCellFromAlias("MaximumFurlAngle"), "=StopFit.Travel")
    # Existing shop calculations express widths as numerical millimetres.
    stop.set(stop.getCellFromAlias("HighEndStopWidth"), "=StopFit.StockWidth / (1 mm)")


if __name__ == "__main__":
    from .load import load_turbine

    _, master = load_turbine(*json.load(sys.stdin))
    dimensions = [master.HighEndStop.MaximumFurlAngle.Value, master.HighEndStop.HighEndStopWidth]
    print("OpenAFPM stop fit: " + json.dumps(dimensions, allow_nan=False), flush=True)
PYTHON_MODULE

cd /app/out

/opt/conda/bin/python - <<'SAVE_DEFAULT'
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile
import FreeCAD as App
from openafpm_cad_core.get_default_parameters import get_default_parameters
from openafpm_cad_core.load import load_turbine

source = Path('/app/inputs/openafpm_cad_core/documents')
destination = Path('openafpm_cad_core/documents').resolve()
config = get_default_parameters('T Shape')
turbine, master = load_turbine(config['magnafpm'], config['furling'], config['user'])

for document in App.listDocuments(True).values():
    if not document.FileName:
        continue
    path = Path(document.FileName)
    relative = path.relative_to(destination)
    if any(not item.isValid() or "Touched" in item.State for item in document.Objects):
        raise RuntimeError(f"Cannot publish an unsettled {document.Label}")
    document.save()
    # A headless document save writes the model and drops the GUI stream.
    # Restore view providers, saved camera and their colors/material resources.
    with ZipFile(source / relative) as original, ZipFile(path, 'a') as saved:
        if 'GuiDocument.xml' in original.namelist():
            gui = original.read('GuiDocument.xml')
            resources = {n.get('file') for n in ElementTree.fromstring(gui).iter() if n.get('file')}
            for name in sorted(resources | {'GuiDocument.xml'}):
                if name not in saved.namelist():
                    saved.writestr(original.getinfo(name), original.read(name))
SAVE_DEFAULT
