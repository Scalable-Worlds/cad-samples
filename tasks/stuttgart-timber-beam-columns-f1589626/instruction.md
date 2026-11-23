Release a glulam layup for a three-member beam-column proof series. One batch of finger-jointed stock must fill all three orders. A, B and C have respective targets (f_req, E_req) of (32, 14200), (30, 13600) and (28, 12600) MPa, using the GL32h, GL30h and GL28h grade designations and mean longitudinal moduli in /app/inputs/Abaqus_timber_beam-column.py.

**Stock and fabrication model**

Treat /app/inputs/data_beam_0.txt as one continuous virtual lamella in file order. Each record is 150 × 120 × 40 mm (length × width × depth). Cut successive 54-record strips and place each intact strip exactly once. A strip may be turned end for end by 180° about Z. As supplied, record numbers increase along +X and the left half lies at lower Y. Across adjoining layers, every pair of finger-joint cell centers must be at least 450 mm apart along X.

Deliver /app/out/layup.vtu as three disjoint conforming linear hexahedral meshes in nominal fabrication coordinates (mm). X runs from 0 to L=8100, +Y across B=120 and +Z down H=480. The section centers (Y,Z) of A, B and C are (0,0), (240,0) and (480,0). Use only VTK_HEXAHEDRON cells, one per 150 × 60 × 40 half-cell, with positive orientation and shared interface nodes within each member. Reference every point and omit ghost and blanking arrays. Integer scalar cell arrays `segment` and `side` give the original `n segment` and its source half (0 left, 1 right). Match vertex coordinates within 0.01 mm per axis. Rate the nominal grid.

The VTU must be a serial UnstructuredGrid with exactly one Piece and pass the `meshio info` command below. `rate.py` must accept the same grid structure in any VTK XML storage, including ASCII, inline binary and appended raw or base64 arrays. Binary data may be uncompressed or compressed with VTK’s zlib, LZ4 or LZMA compressor.

**Proof fixture**

Use N and mm. Pins lie on the support chord at X=0,L. The fixture holds the initial twist and permits bending in both transverse directions. Rate each member as fabricated and after a 180° end-for-end turn about Z. The turn maps an untwisted material point (x,y,z) to (L−x,−y,z), with y,z relative to its section center. Its source segment and side travel with it. The survey stays fixed to the fixture. Use fixture axes for equilibrium. In /app/inputs/data_Imp_meas_0.txt, rows 1 and 3 give lateral bow (m) and twist (rad) at stations iL/100, i=0…100. Rows 2 and 4 give the same quantities halfway between stations. Interpolate linearly through the interleaved samples. Let d(x)=(bow in mm,0). Throughout each installed segment, form r by rotating each point’s installed (y,z) by the twist at that segment’s fixture midpoint, positive from fixture +Y toward +Z.

Use the elastic split-width KaReMo modulus Et_le or Et_ri for the source half in both stress senses. Use the row's ft and fc, or ft_FJ and fc_FJ in a finger-joint cell. The section strain is e(x)−r·u″(x), where u=(u_y,u_z) is additional transverse displacement and e is axial strain at the geometric center. With tension positive, σ is the local modulus times this strain. Under chord compression P and constant applied bending moment m, require at every section
```
∫ σ dA = −P
∫ r σ dA = (0,m) + P (d+u)
```
Use u(0)=u(L)=0 and continuous u and u′. Both integrals cover the full cross-section.

Let W=BH²/6 and I=BH³/12. For each installation at P=68000 N, find f, the largest nonnegative value for which every |m|≤fW keeps −fc≤σ≤ft throughout the member, including both material traces at interfaces. Return f=0 if the preload alone exceeds a strength. From a separate P=0 pure-bending case with m=W × 1 MPa, find E_app=|m|L/(I |u′_z(L)−u′_z(0)|). For each member report the smaller f and the smaller E_app across its two installations.

**Release check**

Run these commands headlessly.
```
/usr/local/bin/meshio info /app/out/layup.vtu
/usr/local/bin/python3 /app/out/rate.py --mesh /app/out/layup.vtu --inputs /app/inputs --out /app/out
```
Save /app/out/rating.json. The program must rate any fabrication-admissible allocation of this stock from the specified mesh and input directory using the preinstalled toolchain. Even for a below-target allocation, it must create the requested output directory if needed, write rating.json there and exit successfully. The JSON object has keys A, B and C, each with exactly the finite numeric fields f_MPa and E_MPa. Allow absolute reporting errors of 0.002 MPa for f and 0.2 MPa for E_app. The submitted mesh’s independently recomputed minimum over all members’ f/f_req and E_app/E_req must be ≥ 1.0009.

Write only your final deliverables under `/app/out`.
