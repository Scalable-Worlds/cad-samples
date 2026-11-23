`/app/inputs/a8.ifc` is a Tekla Structures IFC2X3 export (millimetres, Z up) of a residential block: piled foundation, cast-in-place frame, precast floors, balconies and galleries. Its quantity sets describe each part's own solid, but parts interpenetrate, embed one another and touch, so they cannot be used to order concrete or to schedule surface work. Deliver the take-off program for the cast-in-place concrete.

**Deliverables**
- `/app/out/pour_takeoff.py`, a self-contained script invoked from any directory as `/opt/conda/bin/python /app/out/pour_takeoff.py <model.ifc> <result.json>`. It works from the given model alone, writes the result to the given path, exits with status 0 and uses only the toolchain preinstalled here.
- `/app/out/a8_takeoff.json`: its result for `/app/inputs/a8.ifc`.

**Models.** The script is run on `/app/inputs/a8.ifc` and on other exports of this project that keep its conventions and vocabulary. They may lie anywhere within 1000 km of the origin, at any rotation about the vertical axis, in millimetres or metres, and differ in which parts exist and in how parts are positioned, duplicated or subdivided. In every model poured and void parts have only horizontal and vertical faces (to 0.001°; vertical faces may face any direction); other parts have any shape.

**Parts.** Every solid is a part: an element aggregated by an `IfcElementAssembly`. A part has the first class that applies:
- *Void*: `Name` is `SPARING`. A void is empty space: nothing exists inside it.
- *Annotation*: its material is `MISCELLANEOUS/DUMMY`. Ignore it.
- *Poured*: the assembly's `Cast unit type` (property set `Tekla Assembly`) is `Insitu` and the part is not an `IfcPile`. Its kind is its `Name`, always one of the nine below.
- *Other construction*: every remaining part, piles included.

**Poured body.** The concrete is the union of the poured parts minus other construction and voids. Where kinds overlap, the concrete belongs to the kind that comes first in the pour sequence `POER`, `FUNDERINGSBALK`, `VLOER`, `OPSTORTING`, `KOLOM`, `KOLOMKOP`, `WAND`, `BALK`, `DRUKLAAG`.

**Surface.** Partition the boundary of the poured body at 1 mm resolution: a gap, or a layer of concrete or of other construction, thinner than 1 mm normal to the boundary is absent (so is a stack of such layers), and what lies on either side of it touches. Exports are noisy only below that: apart from slivers of negligible area, no model has such a layer between 0.8 mm and 5 mm thick, and stacked thin layers stay under 0.8 mm in total. Boundary touching poured concrete is an internal joint, not surface. Boundary touching other construction is `contact`. The rest is free surface: `top` facing up, `soffit` facing down, `side` vertical. Each region belongs to the kind that owns the concrete behind it.

**Result.** A JSON object `{"kinds": {"<kind>": {"volume_m3": n, "top_m2": n, "side_m2": n, "soffit_m2": n, "contact_m2": n}}}` with exactly the nine kinds above as keys (zeros where a kind has no concrete).

**Accuracy.** Each volume must be within 0.5 % or 0.05 m³, and each area within 1 % or 0.5 m², whichever is larger, of the exact value for the model as defined here. Curved faces may be faceted with a chordal error up to 1 mm.

Write only your final deliverables under `/app/out`.
