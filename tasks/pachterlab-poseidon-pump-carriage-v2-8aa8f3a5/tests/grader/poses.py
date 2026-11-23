import struct

import numpy as np
from scipy.spatial.transform import Rotation


FACET = np.dtype([
    ("normal", "<f4", (3,)),
    ("vertices", "<f4", (3, 3)),
    ("attribute", "<u2"),
])
SEEDS = {"brace": 7823, "carriage": 7824, "pump-base": 7831, "slide": 7826, "screen-leg": 7832}


def repose(key, source, destination):
    """Rigidly place a whole binary mesh and permute its facets."""
    data = source.read_bytes()
    count = struct.unpack_from("<I", data, 80)[0]
    facets = np.frombuffer(data, FACET, count=count, offset=84).copy()
    rng = np.random.default_rng(SEEDS[key])
    rotation = Rotation.random(random_state=rng).as_matrix()
    translation = rng.uniform(-145, 145, 3)
    facets["vertices"] = facets["vertices"] @ rotation.T + translation
    vertices = facets["vertices"].astype(float)
    normals = np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0])
    facets["normal"] = normals / np.linalg.norm(normals, axis=1)[:, None]
    shuffled = facets[rng.permutation(count)]
    with destination.open("wb") as output:
        output.write(bytes(80))
        output.write(struct.pack("<I", count))
        output.write(shuffled.tobytes())
    return shuffled["vertices"].astype(float), rotation, translation
