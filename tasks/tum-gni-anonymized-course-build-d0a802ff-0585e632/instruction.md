Produce an operating chart for a two-wheel maintenance carrier sharing the lower chord of each roof truss in /app/inputs/model_7_structure.ifc with a suspended service unit. The operator knows the carrier station but may not know where the unit is parked.

Model each IFC2X3 IfcElementAssembly of predefined type TRUSS as a separate vertical plane frame in project-world millimetres, Z up. The placed Axis segments of its aggregated IfcMembers are the centroidal lines. The horizontal chords are continuous through the web joints. Use each placed Body's steel section normal to its axis at analytical midlength, retaining its orientation. Webs are pin-ended axial bars. Use small-displacement linear elasticity, E = 200000 N/mm², Euler–Bernoulli chord bending about the truss-plane normal and continuous lateral restraint. Joints coincide within 0.05 mm.

Measure station along the lower chord from its lower-X end, choosing the lower-Y end if the endpoint X coordinates differ by at most 0.1 mm. Pin this end and use a vertical-reaction roller at the far end. Let L be the analytical span in mm. Apply steel self-weight at 0.0000770 N/mm³ on analytical member lengths, uniformly on each chord and half of each web's weight at each end joint.

The carrier and service unit exert vertical downward point forces on the lower chord axis. At carrier station s the wheels at s and s + 1100 mm carry 65% and 35% of the gross carrier load, in either order. The separate unit applies 40 kN at any station h in [0,L] at least 400 mm from each wheel, including between the wheels wherever both clearances hold. R(s) is the greatest gross carrier load permitted for every such h and both orientations. Include self-weight and the unit simultaneously in every check. If any admissible h violates a limit at zero carrier load, set R(s) to zero.

Limit absolute longitudinal normal stress from axial force and bending throughout the truss to 160 N/mm², and downward displacement anywhere on the lower chord to min(L/600, 20 mm). Limit each web's compressive force to π² E Imin/(2 ℓ²), using its smaller centroidal principal second moment Imin and analytical length ℓ.

Deliver /app/out/rate_trusses.py and /app/out/truss_ratings.json. Run headlessly with

```sh
/opt/conda/bin/python /app/out/rate_trusses.py --structure /app/inputs/model_7_structure.ifc --out /app/out
```

Take both paths from the arguments in any working directory. Create --out and write truss_ratings.json directly there. Reissues contain one or more trusses of the same IFC2X3 family in any vertical plane. The vertically aligned, equal-span chords have a post at every panel station including both ends and one diagonal per bay. Panel runs of 0.7–1.6 m sum to 4–28 m. Chord axes are 1.4–2.4 m apart. Each section uses its supplied chord or web I-profile family scaled uniformly in its plane by a factor from 0.85 to 1.15. Reissued trusses either retain positive capacity throughout or have an admissible permanent overload at every carrier station.

The JSON object has one key "trusses", an array with one record per assembly containing "guid" (assembly GlobalId), "members", "curve" and "witness". In "members" give one record per aggregated IfcMember with "guid" (member GlobalId), "ends_mm" (two world XYZ triples in either order), "A_mm2", "I_mm4" and "c_mm". Report area A. For a chord give centroidal I about the truss-plane normal and c to the farthest in-plane fibre. For a web give Imin and null for c. Allow 1% on A, 2% on I and 0.25 mm on c. Tessellated section curves have at most 0.1 mm chordal deviation.

The curve consists of [s_mm,load_kN] pairs, linear between strictly increasing knots covering 0 through L − 1100. At every station its ordinate must lie between 0.99 R(s) and 1.001 R(s). At a global minimum on that range report a witness with "s_mm", "unit_s_mm" (admissible h), "heavy_at" ("LOW" for the wheel at s, "HIGH" for the other), "mode" ("STRESS", "SAG" or "WEB_BUCKLING"), "member" (member GlobalId) and "xyz_mm" (undeformed world XYZ). Give the governing section's axis point on a chord or the midpoint of a web. For this placement and the charted load, the named demand there must be within 1.5% of its limit. For overload at zero carrier load identify an exceeded limit. World positions and stations have 1 mm tolerance. Any equally governing witness is acceptable.

Write only your final deliverables under `/app/out`.
