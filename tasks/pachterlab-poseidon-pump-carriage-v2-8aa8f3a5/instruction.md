Prepare a mechanical access survey for a side-looking inspection probe on the Poseidon microscope station. The fixture designer has the STL print releases and needs to know which portions of the cylindrical seats the probe window can reach by straight insertion from each end.

Deliver /app/out/bore_access.py. Run it headlessly as `/opt/conda/bin/python /app/out/bore_access.py --mesh MODEL --out DIR`. Save the three surveys below. Each call creates DIR if necessary and writes DIR/survey.json.

| MODEL | DIR |
|---|---|
| /app/inputs/Scope Base v33.stl | /app/out/base |
| /app/inputs/Scope Screen Leg Mount v33.stl | /app/out/screen-leg |
| /app/inputs/Scope Top Stage v33.stl | /app/out/top-stage |

Fresh runs use other whole Poseidon pump and microscope print components. Honor the requested output directory wherever it points and use only MODEL and the preinstalled toolchain. An unreadable MODEL must give a nonzero exit before replacing a report. Each valid input is a single closed, outward-wound binary STL in mm, with a 10–300 mm bounding-box diagonal and coordinates within ±500 mm. File names, triangle order, rotation and translation may vary. The meshes are CAD tessellations. On the cylindrical walls of interest the chordal deviation is at most 0.05 mm and a facet spans at most 16° about the axis.

A land comprises all inward-facing wall on the same nominal circular cylinder, including disconnected portions. Include diameters from 3.5 through 40 mm where the union of occupied azimuths has angular measure of at least 160°. Trace the actual trimmed wall. The same axis can carry several diameters.

Use JSON with a single key `bores`, an unordered array with one object per land, empty when none qualifies. Each object has exactly `c,d,e,r,rows`. In input coordinates, c is the point of the axis closest to the origin, d is a unit vector along the axis in either sense, e is a perpendicular unit vector of your choice for clock zero, and r is the radius in mm. The vectors are XYZ triples. At row k = 0,…,179 use θ = 2k degrees and n = e cos θ + (d × e) sin θ. The wall generator is c + r n + h d, with axial station h in mm.

Each row is [wall, plus, minus]. Each entry is a sorted list of disjoint positive-length [low,high] intervals of h, merging touching intervals and omitting isolated points. Wall gives the stations occupied by the land. Plus and minus give the subsets accessible from the +d and −d ends respectively. The flat-ended probe has a 3.00 mm diameter cylindrical shank extending back to the approach end at infinity. Its axis is parallel to d, at radial vector (r − 1.70 mm)n from the bore axis, leaving 0.20 mm nominal wall clearance. At window station h the tip is 2.00 mm farther in the direction of insertion. Access requires the entire probe to clear the closed solid bounded by the input triangles. Contact blocks insertion. Interval endpoint inclusion has zero weight.

A survey passes when the lands match one-to-one, the Euclidean error in c and the absolute error in r are each at most 0.03 mm, and each unoriented axis is within 0.03°. Require |‖d‖−1|, |‖e‖−1| and |d·e| each ≤0.000001. Reference rows use the nominal cylinder in the reported sense, with e projected onto its normal plane and normalized. For each of the three columns, sum the length of interval symmetric difference over all clock positions and divide by summed reference wall length. Each score on each land must be ≤0.015. All saved surveys and fresh runs must pass.

The following command measures those three scores for any matched land. REF and RUN contain its reference and reported `rows` arrays in the same clocking.

```sh
/opt/conda/bin/python -c '
import json,sys
from shapely import LineString,union_all
def lines(intervals):
    return union_all([LineString([(a,0),(b,0)]) for a,b in intervals])
ref,run=[json.load(open(p)) for p in sys.argv[1:]]
weight=sum(lines(v[0]).length for v in ref)
print([sum((lines(x[j]) ^ lines(y[j])).length for x,y in zip(ref,run,strict=True))/weight for j in range(3)])
' REF RUN
```

Write only your final deliverables under `/app/out`.
