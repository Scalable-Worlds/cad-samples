"""Read the declared PLY records without polygon triangulation or index repair."""

import io
import os
from dataclasses import dataclass, field
from decimal import Decimal

import numpy as np
import trimesh


TYPES = {name: np.dtype(code) for code, names in {
    'i1': ('char', 'int8'), 'u1': ('uchar', 'uint8'),
    'i2': ('short', 'int16'), 'u2': ('ushort', 'uint16'),
    'i4': ('int', 'int32'), 'u4': ('uint', 'uint32'),
    'i8': ('int64',), 'u8': ('uint64',), 'f2': ('float16',),
    'f4': ('float', 'float32'), 'f8': ('double', 'float64'),
}.items() for name in names}
INTEGER_LIMITS = {t: np.iinfo(t) for t in TYPES.values() if t.kind in 'iu'}


@dataclass
class Property:
    name: str
    dtype: np.dtype
    count_type: np.dtype | None = None


@dataclass
class Element:
    name: str
    count: int
    properties: list = field(default_factory=list)


def header_line(stream, ending=None):
    # PLY headers use CR, LF or CRLF. Read exactly to the end of the header,
    # including when the first binary datum happens to be a newline byte.
    line = bytearray()
    while True:
        char = stream.read(1)
        if not char:
            raise ValueError('Missing end of PLY header')
        if char in (b'\r', b'\n'):
            if char == b'\r' and ending != b'\r' and stream.peek(1)[:1] == b'\n':
                char += stream.read(1)
            return line.strip(), char
        line.extend(char)


def read_header(stream):
    signature, ending = header_line(stream)
    if signature != b'ply':
        raise ValueError('Missing PLY signature')
    elements = {}
    encoding = None
    element = None
    while True:
        line, _ = header_line(stream, ending)
        words = line.split()
        if not words or words[0] in (b'comment', b'obj_info'):
            continue
        words = [word.decode('ascii') for word in words]
        if words == ['end_header']:
            break
        if len(words) == 3 and words[0] == 'format' and words[2] == '1.0':
            if encoding is not None:
                raise ValueError('Multiple PLY format declarations')
            encoding = words[1]
        elif len(words) == 3 and words[0] == 'element':
            name, count = words[1], int(words[2])
            if count < 0 or name in elements:
                raise ValueError('Invalid or duplicate PLY element declaration')
            element = elements[name] = Element(name, count)
        elif words[0] == 'property' and element is not None:
            if len(words) == 3:
                prop = Property(words[2], TYPES[words[1]])
            elif len(words) == 5 and words[1] == 'list':
                prop = Property(words[4], TYPES[words[3]], TYPES[words[2]])
                if prop.count_type.kind not in 'iu':
                    raise ValueError('PLY list count must have an integer type')
            else:
                raise ValueError('Malformed PLY property declaration')
            if any(p.name == prop.name for p in element.properties):
                raise ValueError(f'Duplicate {element.name} property {prop.name}')
            element.properties.append(prop)
        else:
            raise ValueError(f'Unexpected PLY header declaration: {words[0]}')
    if encoding not in ('ascii', 'binary_little_endian', 'binary_big_endian'):
        raise ValueError('Missing or unknown PLY encoding')
    return encoding, elements


def ascii_number(token, dtype):
    if dtype.kind in 'iu':
        if not token.lstrip('+-').isascii() or not token.lstrip('+-').isdigit():
            raise ValueError(f'Invalid PLY integer {token}')
        value = int(token)
        bounds = INTEGER_LIMITS[dtype]
        if not bounds.min <= value <= bounds.max:
            raise ValueError(f'Integer outside declared PLY type {dtype}')
        return value
    return dtype.type(token).item()


def ascii_records(stream, element, coordinate_literals):
    for _ in range(element.count):
        line = stream.readline()
        while line and not line.strip():
            line = stream.readline()
        if not line:
            raise ValueError(f'Incomplete PLY {element.name} table')
        fields = iter(line.split())
        row = {}
        literals = {}
        for prop in element.properties:
            try:
                if prop.count_type is None:
                    token = next(fields)
                    value = ascii_number(token, prop.dtype)
                    if element.name == 'vertex' and prop.name in ('x', 'y', 'z'):
                        literals[prop.name] = token if prop.dtype.kind == 'f' else None
                else:
                    count = ascii_number(next(fields), prop.count_type)
                    if count < 0:
                        raise ValueError('Negative PLY list length')
                    value = [ascii_number(next(fields), prop.dtype) for _ in range(count)]
            except StopIteration:
                raise ValueError(f'Truncated PLY {element.name} record') from None
            row[prop.name] = value
        if next(fields, None) is not None:
            raise ValueError(f'Extra fields in PLY {element.name} record')
        if literals:
            coordinate_literals.append(tuple(literals[axis] for axis in 'xyz'))
        yield row


def decimal_precision(literals):
    """Local decimal rounding radii, with absent trailing zeros inferred per field."""
    steps = np.zeros((len(literals), 3))
    for axis in range(3):
        column = [row[axis] for row in literals]
        explicit = [token is not None and ('.' in token or 'e' in token.lower()) for token in column]
        if not any(explicit):
            continue
        values = [Decimal(token) for token in column]
        decimals = [v for v, written in zip(values, explicit) if written]
        finest = min(v.as_tuple().exponent for v in decimals)
        digits = max(len(v.as_tuple().digits) for v in decimals)
        for i, (value, written) in enumerate(zip(values, explicit)):
            exponent = value.as_tuple().exponent
            if not written:
                # Significant-digit formatting omits trailing zeros, including
                # the decimal point on integral results. Infer their resolution
                # from the other literals in this field, never from the units
                # digit of a bare integer alone.
                significant = value.adjusted() + 1 - digits if value else finest
                exponent = min(0, max(finest, significant))
            steps[i, axis] = float(f'5e{exponent - 1}')
    return steps


def binary_records(stream, element, endian, file_size):
    def take(dtype, count):
        size = count * dtype.itemsize
        if count < 0 or size > file_size - stream.tell():
            raise ValueError(f'Truncated or invalid PLY {element.name} record')
        return np.frombuffer(stream.read(size), dtype=dtype.newbyteorder(endian))

    for _ in range(element.count):
        row = {}
        for prop in element.properties:
            if prop.count_type is None:
                value = take(prop.dtype, 1)[0].item()
            else:
                count = int(take(prop.count_type, 1)[0])
                value = take(prop.dtype, count)
            row[prop.name] = value
        yield row


def face_index_name(element):
    properties = {p.name: p for p in element.properties}
    for name in ('vertex_indices', 'vertex_index'):
        if name in properties:
            prop = properties[name]
            break
    else:
        candidates = [p for p in properties.values() if p.count_type is not None and p.dtype.kind in 'iu']
        if len(candidates) != 1:
            raise ValueError('Cannot identify indexed PLY faces')
        prop = candidates[0]
    if prop.count_type is None or prop.dtype.kind not in 'iu':
        raise ValueError('PLY vertex indices must be an integer list')
    return prop.name


def load_indexed_ply(path):
    vertices, faces, coordinate_literals = [], [], []
    coordinate_dtypes = (np.dtype('f8'),) * 3
    with path.open('rb') as stream:
        encoding, elements = read_header(stream)
        ascii_mode = encoding == 'ascii'
        endian = '>' if encoding == 'binary_big_endian' else '<'
        file_size = os.fstat(stream.fileno()).st_size
        data = io.TextIOWrapper(stream, encoding='ascii') if ascii_mode else stream
        for element in elements.values():
            if not element.count:
                continue
            if not element.properties:
                raise ValueError(f'{element.name} records have no declared properties')
            if element.name == 'vertex':
                scalars = {p.name: p.dtype for p in element.properties if p.count_type is None}
                if not all(axis in scalars for axis in 'xyz'):
                    raise ValueError('Missing PLY vertex coordinates')
                coordinate_dtypes = tuple(scalars[axis] for axis in 'xyz')
            if element.name == 'face':
                index_name = face_index_name(element)
            rows = ascii_records(data, element, coordinate_literals) if ascii_mode else binary_records(data, element, endian, file_size)
            for row in rows:
                if element.name == 'vertex':
                    vertices.append(tuple(row[axis] for axis in 'xyz'))
                elif element.name == 'face':
                    indices = row[index_name]
                    if len(indices) != 3:
                        raise ValueError('A PLY face does not declare exactly three indices')
                    faces.append(indices)
        if ascii_mode:
            extra = any(line.strip() for line in data)
        else:
            extra = bool(data.read(1))
        if extra:
            raise ValueError('Data after declared PLY tables')
    coords = np.array(vertices, dtype=float).reshape(-1, 3)
    if not np.isfinite(coords).all():
        raise ValueError('Non-finite PLY vertex coordinates')
    indices = np.array(faces, dtype=np.int64).reshape(-1, 3)
    if np.any(indices < 0) or np.any(indices >= len(coords)):
        raise ValueError('Face index outside the PLY vertex table')
    mesh = trimesh.Trimesh(vertices=coords, faces=indices, process=False)
    mesh.metadata['coordinate_dtypes'] = coordinate_dtypes
    if coordinate_literals:
        mesh.metadata['ascii_coordinate_half_steps'] = decimal_precision(coordinate_literals)
    return mesh
