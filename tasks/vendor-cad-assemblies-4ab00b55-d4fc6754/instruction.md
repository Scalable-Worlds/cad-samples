Prepare replacement laser blanks for a printer frame. The repair shop receives folded sheet-metal parts as STLs and needs a reusable tool that releases the cut contour and signed brake marks together.

Develop on /app/inputs/base.STL, /app/inputs/corner.STL and /app/inputs/left.STL. Their folded CAD masters are /app/inputs/base.STEP, /app/inputs/corner.STEP and /app/inputs/left.STEP. Within each calibration pair the axes agree and the origins may differ by translation.

Each production job supplies a whole frame component or a dimensional variant, in any rigid orientation, as a closed, outward-wound binary STL in millimetres. Every STL vertex coordinate lies in [−500, 500] mm. Coincident facet vertices may differ by up to 0.000001 mm. The bounding-box diagonal is 350–650 mm. The 1.5 mm stock is flat or has planar panels joined by tangent, straight-axis 90° cylindrical bends of 1–2 mm inside radius. The number, size and arrangement of flanges, cutouts and reliefs vary. The sheet develops without overlap into one connected blank, with no two points more than 700 mm apart. Every panel is reachable from the datum panel across at most three bends.

The chord between the axially extreme STL vertices of each panel–bend tangency lies within 0.0005 mm of its registered nominal line even when extended 700 mm at both ends. After registration, chordal deviation is at most 0.09 mm on the stock skins and their trims. The same bound applies to the trim chords after development to mid-gauge. Tessellation preserves the trim loops and the material intervals on each bend's half-angle generator, with interval endpoints within 0.09 mm of nominal along the axis. In the ±0.02 mm strip about a developed half-angle generator, each nominal or tessellated trim branch crosses every parallel generator in the strip once at ≥40° to the axis. The datum is an XYZ point in STL coordinates in the interior of a planar stock face, at least 0.5 mm from its edges. A batch may use the same sheet with datums on different panels and on either side.

Develop the continuous stock skin containing the datum onto the neutral surface halfway through the gauge. Transfer each trim boundary along the local skin normal. Preserve planar dimensions, distances parallel to a bend axis and arc lengths at the neutral radius. Use the neutral-surface image of the datum as the drawing origin. Of the positive STL coordinate axes choose the least aligned with the outward datum normal n, comparing absolute dot products rounded to six decimal places and breaking ties in X, Y, Z order. Over the complete stored datum face, including vertex scatter, the plane minimizing the area integral of squared perpendicular distance has an outward normal within 0.0005° of the nominal normal and selects the same axis. The normalized projection of the chosen axis into the plane gives drawing +X. Drawing +Y is n crossed with drawing +X.

On CUT, draw the blank perimeter and every hole, one simple closed straight-segment LWPOLYLINE per loop, without duplicate cuts. For each bend, mark the generator halfway in angle between the panels, clipped to the material at that angle. Put one LINE per uninterrupted positive-length interval on INNER when the picked skin is concave there and on OUTER when it is convex. Use DXF modelspace XY in mm at 1:1, with zero elevation and zero widths.

Deliver /app/out/cutblank.py and /app/out/blank.dxf. Generate the saved base release with

```sh
/opt/conda/bin/python /app/out/cutblank.py --part /app/inputs/base.STL --datum 161.690406 1.5 156.623074 --out /app/out
```

From any working directory, the same invocation with another STL, datum and output directory must create that directory as needed and write blank.dxf there. Use only the specified STL, datum and preinstalled headless toolchain at runtime.

Reopen the DXF with ezdxf.readfile and compare it in drawing coordinates with the nominal CAD development. For registration, translate the XYZ minima of the FreeCAD Part.read(master).BoundBox to the minima of the STL vertices and project the supplied datum normally onto its corresponding nominal face. The perimeter, holes and signed intervals must match one-to-one, preserving nesting and layer. For each matched pair of curves a,b, require shapely.hausdorff_distance(a,b,densify=0.25) <= 0.15 mm, with nominal curved trim edges sampled at 0.0004 mm chordal deflection. This limit includes source mesh and output chord errors. Every invocation must succeed and every drawing in the batch must meet these conditions.

Write only your final deliverables under `/app/out`.
