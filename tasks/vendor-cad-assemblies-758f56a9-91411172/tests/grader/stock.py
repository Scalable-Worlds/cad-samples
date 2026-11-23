"""Stock geometry, pull coordinates, measurement sites, and nominal steel volume."""

from dataclasses import dataclass
from functools import lru_cache
from math import fsum

import numpy as np
import shapely
import trimesh


@dataclass
class Study:
    basis: np.ndarray
    rectangle_area: float
    volume: float
    lattice: np.ndarray
    cell_points: np.ndarray


@lru_cache(None)
def stock_bytes(source, key):
    """Serialize the same stock for the pull study and independent measurements."""
    if key not in {"lid_placed", "base_placed", "lid_c"}:
        return source.read_bytes()
    mesh = load_stock(source)
    if key == "lid_placed":
        # A quarter turn about Z and a raised datum, exact in binary arithmetic.
        x, y, z = mesh.vertices.T
        mesh.vertices = np.column_stack((-y, x, z + 35.0))
    elif key == "base_placed":
        for angle, axis in ((0.008, (1, 0, 0)), (7, (0, 0, 1))):
            mesh.apply_transform(trimesh.transformations.rotation_matrix(np.radians(angle), axis))
    else:
        # Consistent stock orientation can be inward. Include a cyclic copy
        # of one facet and collinear records supported on a real straight edge.
        mesh.invert()
        triangles = np.asarray(mesh.triangles)
        edges = mesh.vertices[mesh.edges_unique]
        vertical = np.all(edges[:, 0, :2] == edges[:, 1, :2], axis=1)
        a, b = edges[vertical][np.argmax(np.abs(edges[vertical, 1, 2] - edges[vertical, 0, 2]))]
        flat = np.asarray([a, (a + b) / 2, b])
        triangles = np.concatenate((triangles, [triangles[0, [1, 2, 0]], flat,
                                                [a, a, b], [a, a, a]]))
        mesh = trimesh.Trimesh(vertices=triangles.reshape(-1, 3),
                               faces=np.arange(triangles.size // 3).reshape(-1, 3),
                               process=False)
    return mesh.export(file_type="stl")


def load_stock(path):
    """Weld coincident vertices using five-decimal keys, retaining STL coordinates."""
    source = trimesh.load_mesh(path, file_type="stl", process=False)
    coordinates = np.asarray(source.vertices)
    _, first, vertex_id = np.unique(
        np.round(coordinates, 5), axis=0, return_index=True, return_inverse=True
    )
    vertices = coordinates[first]
    faces = vertex_id[source.faces]
    triangles = vertices[faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    faces = faces[np.any(normals != 0, axis=1)]
    _, first_face = np.unique(np.sort(faces, axis=1), axis=0, return_index=True)
    mesh = trimesh.Trimesh(vertices, faces[np.sort(first_face)], process=False)
    mesh.remove_unreferenced_vertices()
    return mesh


def pull_basis(direction):
    """Return the column vectors u, v, d of the published right-handed frame."""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    e = np.eye(3)[np.argmin(np.abs(d))]
    u = np.cross(d, e)
    u /= np.linalg.norm(u)
    return np.column_stack((u, np.cross(d, u), d))


def measurement_sites(polygons, limits):
    start, stop = limits
    pitch = 0.5
    axes = [
        np.arange(lo + offset * pitch, hi, pitch)
        for lo, hi, offset in zip(start, stop, (0.371, 0.613))
    ]
    a, b = np.meshgrid(*axes)
    lattice = np.column_stack((a.ravel(), b.ravel()))

    noded = shapely.union_all(shapely.boundary(polygons), grid_size=1e-8)
    cells = shapely.get_parts(shapely.polygonize(shapely.get_parts(noded)))
    sites = shapely.point_on_surface(cells)
    keep = shapely.distance(sites, shapely.boundary(cells)) > 0.01
    return lattice, shapely.get_coordinates(sites[keep])


def face_heights(triangles, points):
    """Evaluate affine pull depths at centroids of patches inside their triangles."""
    # Near edge-on triangles are retained. Extended precision protects the planar
    # determinant; a patch centroid in a triangle lies within its depth range.
    a, b, c = np.moveaxis(triangles.astype(np.longdouble), 1, 0)
    ba = b - a
    ca = c - a
    pa = points.astype(np.longdouble) - a[:, :2]
    det = ba[:, 0] * ca[:, 1] - ba[:, 1] * ca[:, 0]
    t = np.divide(pa[:, 0] * ca[:, 1] - pa[:, 1] * ca[:, 0], det,
                  out=np.zeros_like(det), where=det != 0)
    s = np.divide(ba[:, 0] * pa[:, 1] - ba[:, 1] * pa[:, 0], det,
                  out=np.zeros_like(det), where=det != 0)
    height = a[:, 2] + t * ba[:, 2] + s * ca[:, 2]
    z = triangles[:, :, 2]
    return np.clip(height.astype(float), z.min(axis=1), z.max(axis=1))


def covered_face_moments(triangles, polygons, exits):
    """Integrate depths of exits with stock above them, or entries with stock below."""
    if len(polygons) < 2:
        return []
    # Two exits along one line imply an intervening entry. The lower exit
    # bounds a trapped interval from below; dually the upper entry bounds it
    # from above. A face can be covered several times, hence the mask union.
    left, right = shapely.STRtree(polygons).query(polygons, predicate="intersects")
    pairs = left < right
    left, right = left[pairs], right[pairs]
    masks = shapely.intersection(polygons[left], polygons[right])
    positive = shapely.area(masks) > 0
    masks, left, right = masks[positive], left[positive], right[positive]
    if not len(masks):
        return []
    centers = shapely.get_coordinates(shapely.centroid(masks))
    hl = face_heights(triangles[left], centers)
    hr = face_heights(triangles[right], centers)
    # Stock facet interiors do not cut through one another, so their order at
    # an overlap centroid holds over that whole convex overlap.
    distinct = hl != hr
    masks = masks[distinct]
    ids = np.where((hl < hr) == exits, left, right)[distinct]
    if not len(masks):
        return []
    order = np.argsort(ids)
    face_ids, first = np.unique(ids[order], return_index=True)
    grouped = np.split(masks[order], first[1:])
    covered = np.array([g[0] if len(g) == 1 else shapely.union_all(g) for g in grouped])
    heights = face_heights(triangles[face_ids], shapely.get_coordinates(shapely.centroid(covered)))
    signed_area = shapely.area(covered) * (-1 if exits else 1)
    return heights * signed_area


def study(mesh, direction):
    """Measure the stated sites and integrate bounded void from the original stock."""
    basis = pull_basis(direction)
    coords = np.asarray(mesh.vertices) @ basis
    limits = np.array((coords[:, :2].min(axis=0), coords[:, :2].max(axis=0)))
    tri = coords[mesh.faces]
    polygons = shapely.polygons(tri[:, :, :2])
    projected = shapely.area(polygons) > 0
    tri, polygons = tri[projected], polygons[projected]
    lattice, cell_points = measurement_sites(polygons, limits)
    # Projected signed area points with d on an outward exit. A whole
    # consistently wound stock may have its orientation reversed.
    exits = shapely.is_ccw(shapely.get_exterior_ring(polygons))
    if mesh.volume < 0:
        exits = ~exits
    moments = []
    for is_exit in (False, True):
        select = exits == is_exit
        moments.extend(covered_face_moments(tri[select], polygons[select], is_exit))
    return Study(basis, float(np.prod(limits[1] - limits[0])), fsum(moments), lattice, cell_points)
