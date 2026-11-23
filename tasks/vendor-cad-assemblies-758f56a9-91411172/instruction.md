We are reviewing the spherical robot's printed shells for production in a two-half mould. For each proposed pull direction, give the toolmaker a solid model of the steel that would be trapped between portions of the part when the halves withdraw in opposite directions.

Use these production STL masters for development.

- /app/inputs/hexapod_spherical_robot/printed_parts/upper_spherical_shell.STL
- /app/inputs/modular_shell/node_shell_upper.STL
- /app/inputs/modular_shell/node_shell_lower.STL

Deliver /app/out/pull.py and the upper spherical panel's study at /app/out/axial.ply, /app/out/side.ply and /app/out/inclined.ply. Generate it with

```sh
/opt/conda/bin/python /app/out/pull.py --part "/app/inputs/hexapod_spherical_robot/printed_parts/upper_spherical_shell.STL" --out /app/out
```

In a fresh box, the same command with another part and destination must write all three files there deterministically using only that STL and the preinstalled headless toolchain. Inputs are spherical skin panels or sensor enclosure halves from this robot, under 200 mm per axis. Their facets, wall profiles, ports, bosses and opening topology can vary. Use millimetres and the original file frame. Treat the stock facets as exact. Weld vertices that coincide to five decimal places in mm and ignore duplicate and zero-area facets. The resulting stock boundary is closed and consistently wound.

The axial, side and inclined setups pull parallel to (0,0,1), (1,0,0) and (2,-1,4) in file XYZ. On each line parallel to the pull, take every bounded air interval between successive intervals of stock. The trapped steel is the closure of the interior of the union of those air intervals. Deliver its boundary as an indexed triangle PLY, allowing disconnected regions. Zero trapped volume requires a PLY with zero faces. Above 0.1 mm³ the PLY must contain faces. At intermediate positive volumes either representation must meet the geometric tolerances.

For each nonempty PLY, require watertightness, consistent outward winding and positive signed volume using the stored vertex indices. Inspect a result with

```sh
/opt/conda/bin/python -c 'import sys,trimesh; m=trimesh.load(sys.argv[1],process=False); print(m.is_watertight,m.is_winding_consistent,m.volume,m.bounds)' /app/out/axial.ply
```

All vertices must be inside the stock's XYZ bounds expanded by 0.05 mm. For unit pull d, use u=unit(d × e) and v=d × u, where e is the first of X,Y,Z minimizing |d·e|. Let A be the area in mm² of the stock's axis-aligned bounding rectangle in (u,v). Nominal volume is the integral of trapped length over that rectangle. Absolute volume error must be at most 0.05 mm × A.

Measure mesh deviation along the pull. On every measured line the total length of the symmetric difference between delivered and nominal intervals must be at most 0.05 mm. Use a square lattice of 0.5 mm pitch inside the projection rectangle, starting (0.371,0.613) pitches from its lower-left corner. Also measure at shapely.point_on_surface of each polygonized cell of the projected stock triangle edges when that point lies more than 0.01 mm from the cell boundary. Use shapely.union_all on the projected boundaries of triangles of nonzero projected area with grid_size=1e-8 mm for polygonization. At degeneracies use solid interiors, merging coincident crossings within 1e-6 mm and discarding tangent contacts. The release scores one if the saved study and all three setups on every rerun pass, and zero otherwise.

Write only your final deliverables under `/app/out`.
