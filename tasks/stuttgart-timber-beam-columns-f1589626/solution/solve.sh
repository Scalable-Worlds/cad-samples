#!/bin/bash
# Regenerate the fabrication issue and rate both fixture installations.
set -euo pipefail
mkdir -p /app/out
rm -f -- /app/out/layup.vtu /app/out/rate.py /app/out/rating.json
work_dir=$(mktemp -d)
trap 'rm -rf -- "$work_dir"' EXIT

cat > "$work_dir/build.py" <<'BUILD_PY'
"""Write the three disjoint members in one VTU Piece. Negative strips turn about Z."""
import csv
from pathlib import Path
import sys

import meshio
import numpy as np


def build(layup, stock_path, destination):
    with stock_path.open(newline="") as stock:
        rows = list(csv.reader(stock))
    columns = [field.strip() for field in rows[0]]
    segment = np.array([int(r[columns.index("n segment")]) for r in rows[1:]],
                       dtype=np.int32).reshape(36, 54)
    joint = np.array([int(r[columns.index("finger joint")]) for r in rows[1:]],
                     dtype=bool).reshape(36, 54)
    assert len(layup) == 3 and all(len(member) == 12 for member in layup)
    assert sorted(abs(strip) for member in layup for strip in member) == list(range(1, 37))

    points = []
    cells = []
    segments = []
    sides = []
    for member_index, member in enumerate(layup):
        # Index the same node from all incident cells of this member.
        vertices = np.arange(len(points), len(points) + 55 * 3 * 13).reshape(55, 3, 13)
        points.extend((x * 150, member_index * 240 + y * 60 - 60, z * 40 - 240)
                      for x, y, z in np.ndindex(vertices.shape))
        previous_joints = []
        for z, strip in enumerate(member):
            source = abs(strip) - 1
            indices = np.arange(54) if strip > 0 else np.arange(53, -1, -1)
            joints = np.flatnonzero(joint[source, indices])
            assert all(abs(a-b) >= 3 for a in joints for b in previous_joints)
            previous_joints = joints
            for x, original in enumerate(indices):
                for y in range(2):
                    cells.append((
                        vertices[x, y, z], vertices[x+1, y, z],
                        vertices[x+1, y+1, z], vertices[x, y+1, z],
                        vertices[x, y, z+1], vertices[x+1, y, z+1],
                        vertices[x+1, y+1, z+1], vertices[x, y+1, z+1],
                    ))
                    segments.append(segment[source, original])
                    sides.append(y if strip > 0 else 1-y)
    meshio.write_points_cells(
        destination, np.array(points, dtype=np.float64),
        [("hexahedron", np.array(cells, dtype=np.int32))],
        cell_data={"segment": [np.array(segments, dtype=np.int32)],
                   "side": [np.array(sides, dtype=np.int32)]},
    )


if __name__ == "__main__":
    # Order is increasing nominal Z for each member, from -240 to +240 mm.
    LAYUP = (
        (-16, 15, -28, -36, 8, -18, -20, -4, 17, 26, -25, 7),
        (-14, -9, -35, 22, -29, 32, -19, -13, 31, -3, 33, 12),
        (-34, 5, 23, -11, -2, 10, -1, -21, -6, 30, 24, -27),
    )
    build(LAYUP, Path("/app/inputs/data_beam_0.txt"), Path(sys.argv[1]))
BUILD_PY

cat > /app/out/rate.py <<'RATE_PY'
#!/usr/local/bin/python3
"""Rate both installations of three pin-ended beam-columns on the nominal grid (N, mm, MPa).

With A_E = integral(E dA), elastic centroid c and centroidal bending
matrix D, eliminating the center strain from the section equations gives

    u'' = -D^-1 [(0,m) + P (d + u + c)]
    sigma = E [-P/A_E - (r-c) . u''].

The two load columns carry preload and the increment for m = W, respectively.
"""
import argparse
import base64
import csv
import json
import lzma
import re
import tempfile
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

import numpy as np

L, B, H = 8100., 120., 480.
DX, DY, DZ = 150., 60., 40.
W, I = B * H**2 / 6, B * H**3 / 12
P = 68000.


def _vtk_dtype(name):
    kind, width = re.fullmatch(r"(Float|Int|UInt)(8|16|32|64)", name).groups()
    return np.dtype({"Float": "f", "Int": "i", "UInt": "u"}[kind] + str(int(width) // 8))


def _vtk_base64(text):
    # VTK writers may encode the header and the payload separately, each
    # with its own base64 padding. Unpadded boundaries need no special case.
    chunks = re.findall(rb"[^=]+={0,2}", b"".join(text.split()))
    return b"".join(base64.b64decode(chunk, validate=True) for chunk in chunks)


def _vtk_lz4(data, size):
    # The numcodecs block format prefixes the LZ4 data with its decoded size.
    from numcodecs.lz4 import decompress
    return decompress(size.to_bytes(4, "little") + data)


def _vtk_payload(data, header_dtype, compressor):
    word = header_dtype.itemsize
    if not compressor:
        length = int(np.frombuffer(data, header_dtype, count=1)[0])
        return data[word:word + length]

    count, size, tail = map(int, np.frombuffer(data, header_dtype, count=3))
    lengths = np.frombuffer(data, header_dtype, count=count, offset=3 * word)
    cursor = (3 + count) * word
    blocks = []
    for i, length in enumerate(lengths):
        block = data[cursor:cursor + int(length)]
        cursor += int(length)
        if compressor == "vtkLZ4DataCompressor":
            block = _vtk_lz4(block, (tail or size) if i == count - 1 else size)
        elif compressor == "vtkZLibDataCompressor":
            block = zlib.decompress(block)
        elif compressor == "vtkLZMADataCompressor":
            block = lzma.decompress(block)
        else:
            raise ValueError(f"Unknown VTU compressor {compressor}")
        blocks.append(block)
    return b"".join(blocks)


def _vtk_document(path):
    """Separate XML metadata from its optional byte-addressed appended stream."""
    from xml.parsers import expat

    xml = Path(path).read_bytes()
    if b"AppendedData" not in xml:
        return ET.fromstring(xml), b"", False

    # A raw payload can contain XML delimiters. Stop at the actual start
    # element, keeping quoted markup, comments and processing instructions
    # distinct from the element that introduces the stream.
    class AppendedSection(Exception):
        pass

    parser = expat.ParserCreate()
    parents = []

    def start_element(name, attributes):
        if len(parents) == 1 and name.split(":")[-1] == "AppendedData":
            raise AppendedSection(parents[0], parser.CurrentByteIndex,
                                  attributes.get("encoding", "base64"))
        parents.append(name)

    parser.StartElementHandler = start_element
    parser.EndElementHandler = lambda name: parents.pop()
    try:
        parser.Parse(xml, True)
    except AppendedSection as section:
        root_name, position, encoding = section.args
        # AppendedData follows the dataset. The prefix contains its full mesh
        # metadata, and the parser never has to consume the binary stream.
        root = ET.fromstring(xml[:position] + f"</{root_name}>".encode("utf-8"))
        if not any(node.get("format") == "appended" for node in root.iter()):
            return root, b"", False
        opening = re.compile(rb"""<(?:[^<>"']|"[^"]*"|'[^']*')*>\s*_""")
        match = opening.match(xml, position)
        if match is None:
            raise ValueError("Missing VTU appended-data marker")
        first = match.end()
        # Encoded offsets count physical bytes (including line breaks).
        # Raw blocks carry their own lengths and may contain a closing tag.
        encoded = encoding == "base64"
        last = xml.index(b"<", first) if encoded else len(xml)
        return root, xml[first:last], encoded
    return ET.fromstring(xml), b"", False


def read_half_cells(path):
    """Return points, hex corner IDs, original segment and source half.

    Read only the arrays the rating uses. Cell and point order are arbitrary.
    Single-piece brick topology puts eight consecutive point IDs in each cell.
    """
    root, appended, encoded = _vtk_document(path)
    piece, = root.findall("{*}UnstructuredGrid/{*}Piece")
    order = ">" if root.get("byte_order") == "BigEndian" else "<"
    header = _vtk_dtype(root.get("header_type", "UInt32")).newbyteorder(order)
    compressor = root.get("compressor")
    offsets = sorted({int(a.get("offset")) for a in root.iter()
                      if a.get("format") == "appended"} | {len(appended)})
    ends = dict(zip(offsets, offsets[1:]))

    def array(location):
        # A serial VTU may store the same fabrication state at several time
        # steps. Read the arrays belonging to the initial state together.
        element = next(a for a in piece.findall(location)
                       if "TimeStep" not in a.attrib or
                       0 in map(int, a.attrib["TimeStep"].split()))
        dtype = _vtk_dtype(element.attrib["type"])
        fmt = element.get("format", "ascii")
        if fmt == "ascii":
            return np.fromstring(element.text, dtype=dtype, sep=" ")
        if fmt == "binary":
            data = _vtk_base64(element.text.encode("ascii"))
        elif fmt == "appended":
            offset = int(element.attrib["offset"])
            data = appended[offset:ends[offset]]
            if encoded:
                data = _vtk_base64(data)
        else:
            raise ValueError(f"Unknown VTU data format {fmt}")
        payload = _vtk_payload(data, header, compressor)
        return np.frombuffer(payload, dtype=dtype.newbyteorder(order))

    points = array("{*}Points/{*}DataArray").reshape(-1, 3)
    cells = array("{*}Cells/{*}DataArray[@Name='connectivity']").reshape(-1, 8)
    segment = array("{*}CellData/{*}DataArray[@Name='segment']").ravel()
    side = array("{*}CellData/{*}DataArray[@Name='side']").ravel()
    return points, cells, segment, side


def read_materials(mesh_path, stock_path):
    with stock_path.open() as file:
        rows = csv.reader(file)
        col = {name.strip(): i for i, name in enumerate(next(rows))}
        stock = {}
        for row in rows:
            suffix = '_FJ' if int(row[col['finger joint']]) else ''
            stock[int(row[col['n segment']])] = tuple(float(row[col[name]]) for name in (
                'Et_le', 'Et_ri', 'ft' + suffix, 'fc' + suffix))

    points, cells, ids, side = read_half_cells(mesh_path)
    # Vertex numbering, including the first local corner, is arbitrary.
    xyz = points[cells].mean(axis=1, dtype=np.float64)
    member = np.rint(xyz[:, 1] / 240).astype(int)
    span = np.rint(xyz[:, 0] / DX - 0.5).astype(int)
    layer = np.rint((xyz[:, 2] + H / 2) / DZ - 0.5).astype(int)
    half = np.rint((xyz[:, 1] - member * 240 + B / 2) / DY - 0.5).astype(int)
    slots = np.ravel_multi_index((member, span, layer, half), (3, 54, 12, 2))
    if len(slots) != 3888 or np.unique(slots).size != len(slots):
        raise ValueError("Expected one cell in every nominal half-cell position")
    values = np.asarray([stock[int(i)] for i in ids])
    grid = np.empty((3, 54, 12, 2, 3))
    grid[member, span, layer, half, 0] = values[np.arange(len(cells)), side]
    grid[member, span, layer, half, 1:] = values[:, 2:]
    return grid


def read_survey(path):
    samples = [np.fromstring(line, sep=',') for line in path.read_text().splitlines()]
    xs = np.linspace(0., L, 201)
    bow, twist = np.empty(201), np.empty(201)
    bow[::2], bow[1::2] = 1000 * samples[0], 1000 * samples[1]
    twist[::2], twist[1::2] = samples[2], samples[3]
    angle = np.interp((np.arange(54) + .5) * DX, xs, twist)
    cuts = np.union1d(xs, np.arange(55) * DX)
    return cuts, np.interp(cuts, xs, bow), angle


def section_properties(modulus, angle):
    """Exact stiffness about each E-weighted centroid, including cell inertias."""
    centers = np.array([[(j + .5) * DY - B / 2, (k + .5) * DZ - H / 2]
                        for k in range(12) for j in range(2)])
    offsets = np.array([[y * DY / 2, z * DZ / 2] for y in (-1, 1) for z in (-1, 1)])
    ec_area = modulus * DY * DZ
    axial = ec_area.sum(axis=1)
    centroid = ec_area @ centers / axial[:, None]
    arm = centers - centroid[:, None, :]
    bending = np.einsum('si,sia,sib->sab', ec_area, arm, arm)
    bending += axial[:, None, None] * np.diag([DY**2 / 12, DZ**2 / 12])
    cs, sn = np.cos(angle), np.sin(angle)
    rotation = np.stack((cs, -sn, sn, cs), axis=1).reshape(-1, 2, 2)
    c = (rotation @ centroid[:, :, None])[:, :, 0]
    stiffness, axes = np.linalg.eigh(rotation @ bending @ rotation.transpose(0, 2, 1))
    # Stress is affine within each half-cell. Keep its four corners, even where
    # another half-cell has the same corner, to evaluate each material trace.
    corner = (arm[:, :, None, :] + offsets).reshape(54, 96, 2)
    principal_arm = corner @ rotation.transpose(0, 2, 1) @ axes
    freq = np.sqrt(P / stiffness)
    stress_map = np.repeat(modulus, 4, axis=1)[:, :, None] * principal_arm * freq[:, None, :]**2
    return axial, c, stiffness, axes, freq, stress_map


def shoot(cuts, bow, span, centroid, axes, frequency):
    """Propagate two load cases plus the two unit initial slopes exactly."""
    length = np.diff(cuts)
    gradient = np.diff(bow) / length
    u = np.zeros((2, 4))
    du = np.zeros((2, 4))
    du[:, 2:] = np.eye(2)
    displacement, velocity = [], []
    for i, s in enumerate(span):
        # For each affine survey piece, the shifted quantity q obeys
        # q'' + P D^-1 q = 0. Restore the physical u,u' before switching pieces.
        shift = np.zeros_like(u)
        slope = np.zeros_like(du)
        shift[:, 0] = centroid[s] + [bow[i], 0.]
        slope[0, 0] = gradient[i]
        shift[1, 1] = W / P
        q = axes[s].T @ (u + shift)
        v = axes[s].T @ (du + slope)
        displacement.append(q)
        velocity.append(v)
        k = frequency[s, :, None]
        cosine, sine = np.cos(k * length[i]), np.sin(k * length[i])
        u = axes[s] @ (cosine * q + sine / k * v) - (shift + length[i] * slope)
        du = axes[s] @ (-k * sine * q + cosine * v) - slope
    # The homogeneous solutions cancel the two right-end displacement residuals.
    weights = np.vstack((np.eye(2), np.linalg.solve(u[:, 2:], -u[:, :2])))
    return np.asarray(displacement) @ weights, np.asarray(velocity) @ weights


def longitudinal_limit(length, k, q, velocity, stress_map, direct, ft, fc):
    """Bracket the full-span capacity, bisecting intervals that can still govern."""
    # Each column of the shifted modal solution is sinusoidal. On an
    # interval of length h, linear interpolation of either stress column
    # has absolute error <= h² max|stress''| / 8. Bound the modal amplitude
    # using |sin(kt)/k| <= t and q'' = -k² q.
    q_bound = np.abs(q) + length[:, None, None] * np.abs(velocity)
    second_bound = np.abs(stress_map) @ (k[:, :, None]**2 * q_bound)

    def min_capacity(reserve, increment):
        # At either moment sign, f*abs(increment) must fit in the
        # smaller of the tension and compression reserves at this trace.
        ratio = np.full_like(reserve, np.inf)
        np.divide(reserve, increment, out=ratio, where=increment > 0.)
        ratio[reserve < 0.] = 0.   # preload alone, even if increment is zero
        return ratio.min(axis=1)

    owner = np.arange(len(length))
    left = np.zeros_like(length)
    right = length.copy()
    upper = np.inf
    while True:
        error = second_bound[owner] * ((right - left)**2 / 8.)[:, None, None]
        lower = np.full(len(owner), np.inf)
        for t in (left, right):
            phase = t[:, None] * k[owner]
            modal = (np.cos(phase)[:, :, None] * q[owner] +
                     (np.sin(phase) / k[owner])[:, :, None] * velocity[owner])
            stress = stress_map[owner] @ modal
            preload = stress[:, :, 0] + direct[owner]
            per_mpa = np.abs(stress[:, :, 1])
            reserve = np.minimum(ft[owner] - preload, fc[owner] + preload)
            upper = min(upper, min_capacity(reserve, per_mpa).min())
            safe = min_capacity(reserve - error[:, :, 0], per_mpa + error[:, :, 1])
            lower = np.minimum(lower, safe)
        # Endpoint stresses bound capacity from above. Subtracting the
        # interpolation allowance for both load columns bounds it from below.
        floor = min(upper, lower.min())
        if upper - floor <= 1.e-9:
            return float(floor)
        active = lower < upper
        mid = (left[active] + right[active]) / 2.
        left, right = np.r_[left[active], mid], np.r_[mid, right[active]]
        owner = np.tile(owner[active], 2)


def rate(material, survey):
    """Rate a single installation, with material indexed in fixture coordinates."""
    cuts, bow, twist = survey
    modulus, ft, fc = material.reshape(54, 24, 3).transpose(2, 0, 1)
    axial, centroid, stiffness, axes, k, stress_map = section_properties(modulus, twist)
    # For P=0, u_z'' = -W (D^-1)_zz. Integrate its exact, constant values.
    compliance_zz = np.sum(axes[:, 1, :]**2 / stiffness, axis=1)
    e_app = L / (I * DX * compliance_zz.sum())
    length = np.diff(cuts)
    span = ((cuts[:-1] + length / 2) / DX).astype(int)
    q, velocity = shoot(cuts, bow, span, centroid, axes, k)
    direct = -P * modulus / axial[:, None]
    direct, ft, fc = [np.repeat(a[span], 4, axis=1) for a in (direct, ft, fc)]
    f = longitudinal_limit(length, k[span], q, velocity, stress_map[span], direct, ft, fc)
    return {'f_MPa': f, 'E_MPa': float(e_app)}


def rate_member(member, survey):
    # Material axes are station, layer, half-width, property. The installed
    # Y halves exchange places, carrying their already assigned moduli and
    # strengths. The fixture survey is applied in its original station order.
    fabricated = rate(member, survey)
    turned = rate(member[::-1, :, ::-1, :], survey)
    return {field: min(value, turned[field]) for field, value in fabricated.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('mesh', 'inputs', 'out'):
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    materials = read_materials(args.mesh, args.inputs / 'data_beam_0.txt')
    survey = read_survey(args.inputs / 'data_Imp_meas_0.txt')
    ratings = {name: rate_member(member, survey) for name, member in zip('ABC', materials)}
    report = json.dumps(ratings, indent=2, allow_nan=False) + '\n'
    args.out.mkdir(parents=True, exist_ok=True)
    # Replace the output entry atomically, including when an older report
    # shares an inode with another file in the release.
    with tempfile.TemporaryDirectory(dir=args.out) as staging:
        pending = Path(staging) / 'rating.json'
        pending.write_text(report, encoding='utf-8')
        pending.replace(args.out / 'rating.json')


if __name__ == '__main__':
    main()
RATE_PY
chmod 755 /app/out/rate.py

/usr/local/bin/python3 "$work_dir/build.py" /app/out/layup.vtu
/usr/local/bin/meshio info /app/out/layup.vtu
/usr/local/bin/python3 /app/out/rate.py --mesh /app/out/layup.vtu --inputs /app/inputs --out /app/out
