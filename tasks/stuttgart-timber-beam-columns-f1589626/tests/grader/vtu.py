"""Read the stored topology and material tuples of a single-piece VTU."""

from base64 import b64encode
import re
from xml.etree import ElementTree as ET
import zlib

import meshio
import numpy as np
from meshio.vtu._vtu import VtuReader, get_grid, numpy_to_vtu_type, vtu_to_numpy_type
from numcodecs import lz4


class ArrayReader(VtuReader):
    """Use the installed byte codecs without converting or filtering cell types."""

    def __init__(self, root):
        # The release CLI treats an empty optional component count as omitted.
        for array in root.findall(".//DataArray[@NumberOfComponents='']"):
            del array.attrib["NumberOfComponents"]
        self.header_type = root.get("header_type", "UInt32")
        self.byte_order = root.get("byte_order")
        self.compression = root.get("compressor")
        self.grid, self.appended_data = get_grid(root)

    def array(self, element):
        if element.get("format", "ascii") == "ascii":
            # Include every token, even values past the final declared offset.
            values = np.asarray((element.text or "").split(),
                                dtype=vtu_to_numpy_type[element.attrib["type"]])
            if "NumberOfComponents" in element.attrib:
                return values.reshape(-1, int(element.attrib["NumberOfComponents"]))
            return values
        return self.read_data(element)


def only(parent, path):
    elements = parent.findall(path)
    assert len(elements) == 1, f"Expected exactly one {path} in {parent.tag}"
    return elements[0]


def scalar(reader, element):
    name = element.get("Name")
    assert int(element.get("NumberOfComponents", 1)) == 1, f"{name} must be scalar"
    values = reader.array(element)
    assert values.dtype.kind in "iu", f"{name} must contain integers"
    return values.ravel()


def material_values(reader, cell_data, name, ncell):
    # VTK TimeStep entries index TimeValues; an unqualified array is static.
    states = set(range(max(1, len(reader.grid.get("TimeValues", "").split()))))
    covered = set()
    reference = None
    for element in cell_data.findall(f"DataArray[@Name='{name}']"):
        selector = element.get("TimeStep", "").split()
        active = states.intersection(map(int, selector)) if selector else states
        if not active:
            continue
        values = scalar(reader, element)
        assert len(values) == ncell, f"Expected one {name} tuple per cell"
        if reference is None:
            reference = values
        else:
            assert np.array_equal(reference, values), f"Conflicting {name} assignments across states"
        covered.update(active)
    assert covered == states, f"Missing {name} assignments at one or more time steps"
    return reference


def xml_root(path):
    try:
        return ET.parse(path).getroot()
    except ET.ParseError:
        # Raw bytes are outside XML syntax. Offsets address the original stream,
        # and several time steps may legally reference the same byte range.
        with open(path, "rb") as file:
            data = file.read()
        opening = re.search(rb"<AppendedData\b[^>]*>\s*_", data)
        assert opening is not None, "Malformed VTU XML"
        end = data.rindex(b"</AppendedData>")
        root = ET.fromstring(data[:opening.end()] + data[end:])
        appended = only(root, "AppendedData")
        assert appended.get("encoding") == "raw", "Malformed encoded VTU stream"
        raw = memoryview(data)[opening.end():end]
        order = ">" if root.get("byte_order") == "BigEndian" else "<"
        word = vtu_to_numpy_type[root.get("header_type", "UInt32")].newbyteorder(order)
        for element in root.findall(".//DataArray[@format='appended']"):
            start = int(element.attrib["offset"])
            assert 0 <= start < len(raw), "Appended array starts outside its data stream"
            first = int(np.frombuffer(raw, word, count=1, offset=start)[0])
            if root.get("compressor"):
                lengths = np.frombuffer(raw, word, count=first, offset=start + 3 * word.itemsize)
                middle = start + (3 + first) * word.itemsize
                stop = middle + sum(map(int, lengths))
                header = b64encode(raw[start:middle])
            else:
                middle = start
                stop = start + word.itemsize + first
                header = b""
            assert stop <= len(raw), "Appended array extends beyond its data stream"
            element.set("format", "binary")
            element.text = (header + b64encode(raw[middle:stop])).decode("ascii")
        root.remove(appended)
        return root


def read_vtu(path):
    """Validate the serialized cell lists before constructing the mesh."""
    root = xml_root(path)
    assert root.tag == "VTKFile" and root.get("type") == "UnstructuredGrid", (
        "Expected a serial VTU UnstructuredGrid"
    )
    reader = ArrayReader(root)
    piece = only(reader.grid, "Piece")
    npoint = int(piece.attrib["NumberOfPoints"])
    ncell = int(piece.attrib["NumberOfCells"])
    points = reader.array(only(only(piece, "Points"), "DataArray"))
    assert points.shape == (npoint, 3), "Points must contain one XYZ tuple per point"

    cells = only(piece, "Cells")
    types, offsets, ids = [scalar(reader, only(cells, f"DataArray[@Name='{name}']"))
                           for name in ("types", "offsets", "connectivity")]
    assert len(types) == len(offsets) == ncell, "Wrong number of cell types or offsets"
    assert np.all(types == 12), "Only VTK_HEXAHEDRON cells may be present"
    assert np.array_equal(offsets, 8 * np.arange(1, ncell + 1)), (
        "Each hexahedron must contain exactly eight connectivity indices"
    )
    assert len(ids) == 8 * ncell, "The final offset must consume the connectivity array"
    assert np.all((ids >= 0) & (ids < npoint)), "Connectivity index outside the point array"
    hexes = ids.astype(np.int64).reshape(-1, 8)
    assert np.all(np.diff(np.sort(hexes, axis=1), axis=1) > 0), "A hexahedron repeats a point ID"
    assert np.unique(hexes).size == npoint, "Every declared point must be referenced by a cell"

    # Hidden points and cells are bits of the reserved vtkGhostType array.
    # vtkGhostLevels is the legacy ghost ownership array.
    for section in ("PointData", "CellData"):
        for attribute in piece.findall(f"{section}/*"):
            assert attribute.get("Name") not in {"vtkGhostType", "vtkGhostLevels"}, (
                f"{section} contains a ghost or blanking array"
            )

    cell_data = only(piece, "CellData")
    labels = {name: [material_values(reader, cell_data, name, ncell)]
              for name in ("segment", "side")}

    for section, count in (("PointData", npoint), ("CellData", ncell)):
        for element in piece.findall(f"{section}/DataArray"):
            name = element.get("Name")
            if section == "CellData" and name in labels:
                continue
            assert len(reader.array(element)) == count, f"Wrong {section} tuple count for {name}"

    return meshio.Mesh(points, [("hexahedron", hexes)], cell_data=labels)


def write_appended(path, mesh, encoding, *, compressor="zlib"):
    """Write a single VTK piece with byte offsets into a raw or base64 data stream."""
    assert encoding in ("raw", "base64")
    codec = {"zlib": "vtkZLibDataCompressor", "lz4": "vtkLZ4DataCompressor"}[compressor]
    root = ET.Element("VTKFile", type="UnstructuredGrid", version="1.0",
                      byte_order="LittleEndian", header_type="UInt32", compressor=codec)
    hexes = mesh.cells[0].data
    grid = ET.SubElement(root, "UnstructuredGrid")
    piece = ET.SubElement(grid, "Piece", NumberOfPoints=str(len(mesh.points)),
                          NumberOfCells=str(len(hexes)))
    stream = bytearray()

    def append_array(section, name, values):
        values = np.asarray(values)
        vtk_type = numpy_to_vtu_type[values.dtype.newbyteorder("=")]
        element = ET.SubElement(section, "DataArray", Name=name, format="appended",
                                type=vtk_type, offset=str(len(stream)))
        if values.ndim == 2:
            element.set("NumberOfComponents", str(values.shape[1]))
        data = values.astype(values.dtype.newbyteorder("<"), copy=False).tobytes()
        block_size = 32768
        blocks = [data[start:start + block_size] for start in range(0, len(data), block_size)]
        # VTK stores bare LZ4 blocks; numcodecs prefixes each with a four-byte size.
        compressed = ([zlib.compress(block) for block in blocks] if compressor == "zlib" else
                      [lz4.compress(block)[4:] for block in blocks])
        header = np.asarray([len(blocks), block_size, len(data) % block_size,
                             *map(len, compressed)], dtype="<u4").tobytes()
        payload = b"".join(compressed)
        for part in (header, payload):
            stream.extend(part if encoding == "raw" else b64encode(part))

    append_array(ET.SubElement(piece, "Points"), "Points", mesh.points)
    cells = ET.SubElement(piece, "Cells")
    append_array(cells, "connectivity", hexes.ravel())
    append_array(cells, "offsets", np.arange(1, len(hexes) + 1, dtype=np.int64) * 8)
    append_array(cells, "types", np.full(len(hexes), 12, dtype=np.uint8))
    cell_data = ET.SubElement(piece, "CellData")
    for name in ("segment", "side"):
        append_array(cell_data, name, mesh.cell_data[name][0])
    ET.SubElement(root, "AppendedData", encoding=encoding).text = "_"
    prefix, suffix = ET.tostring(root).split(b"</AppendedData>")
    with open(path, "wb") as file:
        file.write(prefix)
        file.write(stream)
        file.write(b"\n</AppendedData>")
        file.write(suffix)
