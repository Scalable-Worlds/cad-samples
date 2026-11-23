from contextlib import contextmanager, ExitStack
from decimal import Decimal
from itertools import permutations, product
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

import meshio
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grader.allocations import challenge_layout
from grader.mechanics import (
    RELEASE_RATIO, meets_target, rate_member, read_stock, read_survey, within_report,
)
from grader.vtu import read_vtu, write_appended


MEMBERS = ("A", "B", "C")
TARGETS = {"A": (32.0, 14200.0), "B": (30.0, 13600.0), "C": (28.0, 12600.0)}
OUTPUT = Path("/app/out")
PUBLIC_INPUTS = Path("/app/inputs")
SUPERVISOR = Path(__file__).parent / "grader" / "supervise.py"


@pytest.fixture(scope="module")
def stock():
    return read_stock(PUBLIC_INPUTS)


def load_mesh(path):
    assert path.is_file(), "Missing fabrication mesh layup.vtu"
    try:
        return read_vtu(path)
    except (Exception, SystemExit) as exc:
        assert False, f"layup.vtu cannot be read as a VTU mesh: {exc}"


@pytest.fixture(scope="module")
def mesh():
    return load_mesh(OUTPUT / "layup.vtu")


def test_linear_hexahedra(mesh):
    """Use one VTK hexahedron per 150 × 60 × 40 half-cell in all three members."""
    assert len(mesh.cells[0].data) == 3 * 54 * 2 * 12, "Wrong number of fabrication half-cells"


def nominal_lattice(mesh):
    xyz = mesh.points
    assert np.all(np.isfinite(xyz)), "Half-cell nodes must have finite X,Y,Z positions"
    axes = (150.0 * np.arange(55),
            np.array([-60., 0., 60., 180., 240., 300., 420., 480., 540.]),
            40.0 * np.arange(13) - 240.)
    indices = []
    for i, axis in enumerate(axes):
        ix = np.abs(xyz[:, i, None] - axis).argmin(axis=1)
        target = axis[ix]
        assert np.all((xyz[:, i] >= target - 0.01) & (xyz[:, i] <= target + 0.01)), (
            f"Node coordinate along {'XYZ'[i]} is outside its nominal ±0.01 mm interval"
        )
        indices.append(ix)
    grid = np.column_stack(indices)
    return grid, grid[mesh.cells[0].data]


@pytest.fixture(scope="module")
def lattice(mesh):
    return nominal_lattice(mesh)


def test_nominal_member_grids(lattice):
    """Match the three nominal fabrication grids within 0.01 mm per axis."""
    nodes, _ = lattice
    found = np.unique(nodes, axis=0)
    assert len(found) == 3 * 55 * 3 * 13, "Some nominal grid vertices are absent"


def test_welded_interfaces(lattice):
    """Use shared interface nodes within each conforming member."""
    nodes, _ = lattice
    assert len(np.unique(nodes, axis=0)) == len(nodes), (
        "Two distinct connected node IDs occupy the same nominal grid vertex"
    )


def brick_bounds(lattice):
    _, corners = lattice
    lower = corners.min(axis=1)
    offset = corners - lower[:, None, :]
    assert np.all(offset <= 1), "A hexahedron spans more than one nominal grid interval"
    assert np.all(lower[:, 1] % 3 < 2), "A cell bridges the gap between disjoint members"
    masks = offset @ np.array([1, 2, 4])
    assert np.all(np.sort(masks, axis=1) == np.arange(8)), (
        "Each hexahedron must cover all eight corners of one half-cell"
    )
    return lower, offset


@pytest.fixture(scope="module")
def boxes(lattice):
    return brick_bounds(lattice)


def test_conforming_cell_coverage(boxes):
    """Fill every half-cell of each disjoint member exactly once."""
    lower, _ = boxes
    assert len(lower) == len(np.unique(lower, axis=0)) == 3 * 54 * 2 * 12, (
        "Duplicated or uncovered material half-cells"
    )


def check_orientation(offset):
    # VTK local edges from vertex 0 end at 1, 3, 4. Every orientation of
    # this rectangular brick is covered by the signed axis permutations.
    # At ±0.01 mm the vertex perturbations cannot reverse the Jacobian of
    # an oriented brick whose shortest nominal edge is 40 mm.
    frame = offset[:, [1, 3, 4]] - offset[:, :1]
    unit = np.array([[0,0,0], [1,0,0], [1,1,0], [0,1,0],
                     [0,0,1], [1,0,1], [1,1,1], [0,1,1]])
    assert np.array_equal(unit @ frame + offset[:, :1], offset), (
        "VTK faces and edges do not bound their nominal brick"
    )
    assert np.all(np.linalg.det(frame) > 0), "Inverted hexahedral node order"


def test_oriented_hexahedral_faces(boxes):
    """Use positively oriented linear hexahedra on each nominal half-cell."""
    check_orientation(boxes[1])


def stock_halves(mesh, stock):
    supplied = mesh.cell_data["segment"][0]
    sides = mesh.cell_data["side"][0]
    row_of = {int(segment): i for i, segment in enumerate(stock["segment"])}
    assert all(int(s) in row_of for s in supplied), "Unknown source segment"
    assert np.all((sides == 0) | (sides == 1)), "Source side must be 0 or 1"
    rows = np.fromiter((row_of[int(s)] for s in supplied), dtype=np.int64)
    return np.column_stack((rows, sides.astype(np.int64)))


@pytest.fixture(scope="module")
def cell_sources(mesh, stock):
    return stock_halves(mesh, stock)


def check_custody(cell_sources, stock):
    source = 2 * cell_sources[:, 0] + cell_sources[:, 1]
    assert np.array_equal(np.sort(source), np.arange(2 * len(stock["segment"]))), (
        "Missing or duplicated source half in the fabrication batch"
    )


def test_batch_stock_custody(cell_sources, stock):
    """Place each strip exactly once, retaining both original source halves."""
    check_custody(cell_sources, stock)


def material_layout(lower, cell_sources):
    assert len(lower) == len(np.unique(lower, axis=0)) == 3 * 54 * 2 * 12, (
        "Cannot assign materials without one hexahedron per nominal half-cell"
    )
    x, y, z = lower.T
    layout = np.empty((3, 54, 2, 12, 2), dtype=np.int64)
    layout[y // 3, x, y % 3, z] = cell_sources
    return layout


@pytest.fixture(scope="module")
def allocation(boxes, cell_sources):
    return material_layout(boxes[0], cell_sources)


def check_strips(allocation, member):
    rows = allocation[member, ..., 0]
    assert np.array_equal(rows[:, 0], rows[:, 1]), "Two half-widths came from different records"
    sequence = rows[:, 0]
    start = (sequence[0] // 54) * 54
    forward = start + np.arange(54)[:, None]
    assert np.all(np.all(sequence == forward, axis=0) | np.all(sequence == forward[::-1], axis=0)), (
        "A layer is not one complete strip in forward or reversed file order"
    )


def check_sides(allocation, member):
    rows = allocation[member, :, 0, :, 0]
    flipped = rows[0] > rows[-1]
    sides = allocation[member, ..., 1]
    assert np.all(sides[:, 0] == flipped) and np.all(sides[:, 1] == ~flipped), (
        "Source side disagrees with the orientation of its intact strip"
    )


def check_stagger(allocation, stock, member):
    rows = allocation[member, :, 0, :, 0]
    joint = np.asarray(stock["joint"])[rows]
    locations = [np.flatnonzero(joint[:, layer]) for layer in range(12)]
    for layer, (upper, lower) in enumerate(zip(locations, locations[1:])):
        separations = 150 * np.abs(upper[:, None] - lower)
        assert np.all(separations >= 450), f"Finger joints too close in layers {layer+1} and {layer+2}"


@pytest.mark.parametrize("member", range(3), ids=list("ABC"))
def test_intact_strips(allocation, member):
    """Cut successive groups of 54 records into intact full-width strips."""
    check_strips(allocation, member)


@pytest.mark.parametrize("member", range(3), ids=list("ABC"))
def test_end_turn_exchanges_source_sides(allocation, member):
    """A 180° turn about Z reverses record order and swaps left and right."""
    check_sides(allocation, member)


@pytest.mark.parametrize("member", range(3), ids=list("ABC"))
def test_adjoining_finger_joints_are_staggered(allocation, stock, member):
    """Keep every directly adjoining finger-joint pair at least 450 mm apart along X."""
    check_stagger(allocation, stock, member)


def distinct_object(pairs):
    record = {}
    for key, value in pairs:
        if key in record:
            raise ValueError(f"Duplicate JSON key {key!r}")
        record[key] = value
    return record


def invalid_constant(token):
    raise ValueError(f"{token} is not a JSON number")


def decode_report(data):
    try:
        # Let the JSON decoder recognize the permitted Unicode encodings.
        report = json.loads(data, object_pairs_hook=distinct_object,
                            parse_float=Decimal, parse_int=Decimal, parse_constant=invalid_constant)
    except (ValueError, OSError, UnicodeError) as exc:
        assert False, f"Cannot read rating.json: {exc}"
    return report


@contextmanager
def generated_file(requested, name, boundary):
    """Open a staged output through its contained, regular-file target."""
    try:
        directory = requested.resolve(strict=True)
        path = (directory / name).resolve(strict=True)
        assert directory.is_relative_to(boundary) and path.is_relative_to(directory), (
            f"Generated {name} must resolve inside the requested output directory"
        )
        # Pin the resolved directories while opening. A changing link or a
        # special file must never redirect or block the privileged reader.
        parts = path.relative_to(boundary).parts
        assert parts, f"{name} must name a file"
        with ExitStack() as stack:
            flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
            parent = os.open(boundary, flags)
            stack.callback(os.close, parent)
            for part in parts[:-1]:
                parent = os.open(part, flags, dir_fd=parent)
                stack.callback(os.close, parent)
            descriptor = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=parent)
            stack.callback(os.close, descriptor)
            assert stat.S_ISREG(os.fstat(descriptor).st_mode), f"Generated {name} must name a regular file"
            yield Path(f"/proc/self/fd/{descriptor}")
    except (OSError, RuntimeError) as exc:
        assert False, f"Rating command did not publish a readable {name}: {exc}"


def assert_report_schema(report):
    assert isinstance(report, dict) and set(MEMBERS) <= report.keys(), "Expected ratings keyed by A, B and C"
    for name in MEMBERS:
        record = report[name]
        assert isinstance(record, dict) and set(record) == {"f_MPa", "E_MPa"}, f"{name}: invalid rating fields"
        for field, value in record.items():
            assert isinstance(value, Decimal) and value.is_finite(), f"{name}.{field} must be a finite number"


@pytest.fixture(scope="module")
def survey():
    return read_survey(PUBLIC_INPUTS)


@pytest.fixture(scope="module")
def issued_report():
    path = OUTPUT / "rating.json"
    assert path.is_file(), "Missing issued rating.json"
    try:
        return decode_report(path.read_bytes())
    except OSError as exc:
        assert False, f"Cannot read issued rating.json: {exc}"


@pytest.fixture(scope="module")
def issued_values(issued_report):
    assert_report_schema(issued_report)
    return issued_report


@pytest.fixture(scope="module")
def issued_ratings(allocation, stock, survey):
    return {name: rate_member(member, stock, survey) for name, member in zip(MEMBERS, allocation)}


def test_rating_report_structure(issued_report):
    """“The JSON object has keys A, B and C, each containing exactly the finite numeric fields f_MPa and E_MPa.”"""
    assert_report_schema(issued_report)


@pytest.mark.parametrize("name", MEMBERS)
def test_issued_bending_capacity(name, issued_values, issued_ratings):
    """Report the smaller proof strength of the two installations at the specified proof preload."""
    assert within_report(issued_ratings[name], 0, issued_values[name]["f_MPa"]), f"Wrong proof strength for {name}"


@pytest.mark.parametrize("name", MEMBERS)
def test_issued_pure_bending_modulus(name, issued_values, issued_ratings):
    """Report the smaller apparent modulus from the two installed pure-bending cases."""
    assert within_report(issued_ratings[name], 1, issued_values[name]["E_MPa"]), f"Wrong apparent modulus for {name}"


@pytest.mark.parametrize("name", MEMBERS)
def test_member_strength_qualifies_for_order(name, issued_ratings):
    """Each member's independently recomputed f/f_req must reach the release ratio."""
    actual = issued_ratings[name][0]
    target = TARGETS[name][0]
    assert meets_target(issued_ratings[name], 0, target), (
        f"{name} proof-strength ratio {actual / target:.9f} falls below {RELEASE_RATIO} "
        f"({actual:.7f} MPa against a grade target of {target:g} MPa)"
    )


@pytest.mark.parametrize("name", MEMBERS)
def test_member_stiffness_qualifies_for_order(name, issued_ratings):
    """Each member's independently recomputed E_app/E_req must reach the release ratio."""
    actual = issued_ratings[name][1]
    target = TARGETS[name][1]
    assert meets_target(issued_ratings[name], 1, target), (
        f"{name} apparent-modulus ratio {actual / target:.9f} falls below {RELEASE_RATIO} "
        f"({actual:.7f} MPa against a grade target of {target:g} MPa)"
    )


def stage_release(source):
    # Copy independent inodes, retaining internal hard links and user/group ownership.
    subprocess.run(["cp", "-a", "--", str(source), str(OUTPUT)], check=True)
    for path in (OUTPUT, *OUTPUT.rglob("*")):
        info = path.lstat()
        if info.st_uid == info.st_gid == 65534 and not stat.S_ISLNK(info.st_mode):
            # The entry runner normalizes ownership before pytest starts. A helper
            # formerly executable through group/other may now use its owner bits.
            mode = stat.S_IMODE(info.st_mode)
            access = ((mode & 0o050) << 3) | ((mode & 0o005) << 6)
            if path == OUTPUT:
                access |= 0o700
            if (mode | access) != mode:
                path.chmod(mode | access)


@contextmanager
def working_release(resume=None):
    """Run at the published path while the original tree is privately held intact."""
    # Both directories live under the same protected parent. Moving the original
    # keeps every inode, mode, ownership and hard-link relationship for restoration.
    with tempfile.TemporaryDirectory(dir=OUTPUT.parent) as private:
        original = Path(private) / "original"
        OUTPUT.rename(original)
        try:
            if resume is None:
                stage_release(original)
            else:
                resume.rename(OUTPUT)
            yield
        finally:
            if OUTPUT.exists():
                OUTPUT.rename(Path(private) / "finished")
            original.rename(OUTPUT)


@contextmanager
def scratch_directory():
    with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
        path = Path(tmp)
        os.chown(path, 65534, 65534)
        yield path


def run_delivered(command, scratch):
    env = {**os.environ, "HOME": str(scratch), "TMPDIR": str(scratch), "MPLCONFIGDIR": str(scratch)}
    env.pop("PYTEST_CURRENT_TEST", None)
    with tempfile.TemporaryFile() as out, tempfile.TemporaryFile() as err:
        return subprocess.run([sys.executable, "-I", "-S", str(SUPERVISOR), *command],
                              cwd=scratch, env=env, stdin=subprocess.DEVNULL,
                              stdout=out, stderr=err).returncode


def test_headless_mesh_info():
    """The issued layup must complete the published meshio info release command."""
    path = OUTPUT / "layup.vtu"
    assert path.is_file(), "Missing fabrication mesh layup.vtu"
    with working_release(), scratch_directory() as scratch:
        status = run_delivered(["/usr/local/bin/meshio", "info", str(path)], scratch)
    assert status == 0, f"meshio info exited with status {status} on the issued VTU"


def brick_orientations(corners):
    """Enumerate the orientation-preserving permutations of VTK brick corners."""
    lookup = {tuple(point): i for i, point in enumerate(corners)}
    orders = []
    for axes in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            frame = np.eye(3, dtype=int)[list(axes)] * np.array(signs)[:, None]
            if np.linalg.det(frame) > 0:
                xyz = corners @ frame
                xyz -= xyz.min(axis=0)
                orders.append([lookup[tuple(p)] for p in xyz])
    return np.asarray(orders)


def write_challenge(path, layout, stock, case):
    """Write nominal hexahedra carrying the input record ID and original half."""
    points = []
    hexes = []
    labels = []
    offsets = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    cell_corners = np.indices((54, 2, 12)).reshape(3, -1).T
    grid_dims = (55, 3, 13)
    grid_nodes = np.indices(grid_dims).reshape(3, -1).T
    for index, member in enumerate(layout):
        coords = grid_nodes * [150.0, 60.0, 40.0] + [0.0, index * 240.0 - 60.0, -240.0]
        corners = cell_corners[:, None, :] + offsets
        connectivity = np.ravel_multi_index(corners.transpose(2, 0, 1), grid_dims)
        hexes.append(len(grid_nodes) * index + connectivity)
        points.append(coords)
        labels.append(member.reshape(-1, 2))
    points, hexes, labels = map(np.concatenate, (points, hexes, labels))
    rng = np.random.default_rng(1)
    rotations = brick_orientations(offsets)
    rotation = rng.permutation(len(hexes)) % len(rotations)
    hexes = np.take_along_axis(hexes, rotations[rotation], axis=1)
    # Move both external faces along each physical axis. A contraction and
    # an expansion encode the same nominal grid and point/cell incidence.
    face_offset = (-0.009, 0.009)[case % 2]
    for axis, relative, extent in (
            (0, points[:, 0] - 4050., 4050.),
            (1, points[:, 1] - 240. * np.rint(points[:, 1] / 240.), 60.),
            (2, points[:, 2], 240.)):
        exterior = np.abs(relative) == extent
        points[exterior, axis] += face_offset * np.sign(relative[exterior])
    middle = ((grid_nodes[:, 0] > 0) & (grid_nodes[:, 0] < 54)
              & (grid_nodes[:, 1] == 1) & (grid_nodes[:, 2] > 0) & (grid_nodes[:, 2] < 12))
    middle = np.tile(middle, len(layout))
    points[middle] += rng.choice([-0.005, 0.005], size=(middle.sum(), 3))
    order = rng.permutation(len(points))
    cells = rng.permutation(len(hexes))
    side = labels[cells, 1]
    segment = stock["segment"][labels[cells, 0]]
    if case % 2 == 0:
        # Two-dimensional scalar columns explicitly declare one VTK component.
        side = side[:, None]
        segment = segment[:, None]
    mesh = meshio.Mesh(points[order], [("hexahedron", np.argsort(order)[hexes[cells]])],
                       cell_data={"side": [side], "segment": [segment]})
    if case < 2:
        mesh.write(path, binary=bool(case))
    else:
        write_appended(path, mesh, ("raw", "base64")[case - 2])


@contextmanager
def input_directory(scratch, relocate):
    if not relocate:
        yield PUBLIC_INPUTS
        return
    selected = scratch / "inputs"
    shutil.copytree(PUBLIC_INPUTS, selected)
    original_mode = stat.S_IMODE(PUBLIC_INPUTS.stat().st_mode)
    try:
        PUBLIC_INPUTS.chmod(0o700)
        yield selected
    finally:
        PUBLIC_INPUTS.chmod(original_mode)


def run_rating(mesh_path, inputs, requested, scratch):
    status = run_delivered(["/usr/local/bin/python3", str(OUTPUT / "rate.py"),
                            "--mesh", str(mesh_path), "--inputs", str(inputs),
                            "--out", str(requested)], scratch)
    assert status == 0, f"Rating command exited with status {status}"
    with generated_file(requested, "rating.json", OUTPUT if requested == OUTPUT else scratch) as path:
        report = decode_report(path.read_bytes())
    assert_report_schema(report)
    return report


CHALLENGES = ("restacked", "reoriented", "ordinary-tension", "interface-trace")


@pytest.fixture(scope="module")
def challenge_cases(stock, survey):
    cases = []
    for case in range(len(CHALLENGES)):
        layout = challenge_layout(stock["joint"], case)
        ratings = {name: rate_member(member, stock, survey) for name, member in zip(MEMBERS, layout)}
        cases.append((layout, ratings))
    return cases


@pytest.fixture(scope="module")
def issued_replay():
    # Retain the actual executed installation, even if the program legitimately
    # reorganizes dependencies. The existing-directory replay resumes this tree.
    for name in ("rate.py", "layup.vtu", "rating.json"):
        assert (OUTPUT / name).is_file(), f"Missing issued {name}"
    with tempfile.TemporaryDirectory(dir=OUTPUT.parent) as private:
        executed = Path(private) / "executed"
        with working_release():
            with scratch_directory() as scratch:
                report = run_rating(OUTPUT / "layup.vtu", PUBLIC_INPUTS, OUTPUT, scratch)
            with generated_file(OUTPUT, "layup.vtu", OUTPUT) as path:
                mesh = load_mesh(path)
            with scratch_directory() as scratch:
                info = run_delivered(["/usr/local/bin/meshio", "info", str(OUTPUT / "layup.vtu")], scratch)
            OUTPUT.rename(executed)
        yield mesh, report, info, executed


def test_final_issued_mesh_info(issued_replay):
    """The mesh left by the release commands must still pass headless meshio info."""
    assert issued_replay[2] == 0, f"meshio info exited with status {issued_replay[2]} on the final issued VTU"


def test_issued_command_preserves_fabrication(issued_replay, allocation, stock):
    """The executed release still delivers the conforming mesh with its intact, staggered stock."""
    mesh = issued_replay[0]
    lattice = nominal_lattice(mesh)
    nodes = lattice[0]
    assert len(nodes) == len(np.unique(nodes, axis=0)) == 3 * 55 * 3 * 13, (
        "Regenerated mesh has missing or unwelded fabrication vertices"
    )
    lower, offset = brick_bounds(lattice)
    check_orientation(offset)
    sources = stock_halves(mesh, stock)
    check_custody(sources, stock)
    final = material_layout(lower, sources)
    for member in range(3):
        check_strips(final, member)
        check_sides(final, member)
        check_stagger(final, stock, member)
    assert np.array_equal(final, allocation), "Release command changed the supplied material assignment"


@pytest.fixture(scope="module", params=("issued", *CHALLENGES))
def rated_case(request, stock, challenge_cases):
    if request.param == "issued":
        return request.getfixturevalue("issued_replay")[1], request.getfixturevalue("issued_ratings")
    case = CHALLENGES.index(request.param)
    layout, expected = challenge_cases[case]
    # Continue from the issued run with its complete, valid report still in
    # place. Every other allocation starts with an independent installation.
    executed = request.getfixturevalue("issued_replay")[3] if case == 0 else None
    assert (OUTPUT / "rate.py").is_file(), "Missing issued rating program"
    with working_release(resume=executed), scratch_directory() as scratch:
        mesh_path = scratch / "stock.vtu"
        write_challenge(mesh_path, layout, stock, case)
        requested = OUTPUT if case == 0 else scratch / "ratings"
        with input_directory(scratch, relocate=(case == 1)) as inputs:
            report = run_rating(mesh_path, inputs, requested, scratch)
        return report, expected


@pytest.mark.parametrize("name", MEMBERS)
def test_rating_command_proof_capacity(rated_case, name):
    """Rate each admissible allocation for both installation orientations, including below-target stock maps."""
    report, expected = rated_case
    assert within_report(expected[name], 0, report[name]["f_MPa"]), f"Command returned wrong proof strength for {name}"


@pytest.mark.parametrize("name", MEMBERS)
def test_rating_command_effective_modulus(rated_case, name):
    """Return the lower installed pure-bending modulus independently of which installation limits strength."""
    report, expected = rated_case
    assert within_report(expected[name], 1, report[name]["E_MPa"]), f"Command returned wrong apparent modulus for {name}"


def test_rating_command_reads_appended_lz4(mesh, issued_ratings):
    """Accept the same issued grid with VTK LZ4-compressed, base64-appended arrays."""
    with working_release(), scratch_directory() as scratch:
        mesh_path = scratch / "stock.vtu"
        write_appended(mesh_path, mesh, "base64", compressor="lz4")
        report = run_rating(mesh_path, PUBLIC_INPUTS, scratch / "ratings", scratch)
    for name, expected in issued_ratings.items():
        assert within_report(expected, 0, report[name]["f_MPa"]), f"Wrong strength for {name} on the LZ4 grid"
        assert within_report(expected, 1, report[name]["E_MPa"]), f"Wrong modulus for {name} on the LZ4 grid"
