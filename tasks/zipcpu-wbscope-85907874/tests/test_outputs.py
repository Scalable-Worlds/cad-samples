import atexit
import ctypes
import glob
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from functools import cache

import pytest


OUT = Path("/app/out")
INPUTS = Path("/app/inputs")
GRADER = Path(__file__).parent / "grader"
LIBRARY = ("memscope.v", "memscopc.v", "axisrle.v", "skidbuffer.v", "sfifo.v")
BUILDS = ((16, 32, 5), (16, 64, 5), (32, 128, 9))
DOMAIN_BUILDS = BUILDS + ((12, 128, 9), (12, 32, 2))
WIDTH_SCENARIOS = (2, 3, 4, 6, 9, 15, 21, 27)


def build_id(build):
    a, d, _ = build
    return f"addr{a}-data{d}"


UNPRIVILEGED = dict(user=65534, group=65534, extra_groups=[])


@pytest.fixture(scope="module")
def sources():
    """Preserve the complete internal include tree when relocating the library."""
    rtl = OUT / "rtl"
    assert rtl.is_dir(), f"Missing RTL directory {rtl}"
    files = {}
    def copy_directory(path, ancestors):
        resolved = path.resolve(strict=True)
        assert resolved.is_relative_to(rtl), f"External library directory: {path}"
        assert resolved not in ancestors, f"Cyclic include directory: {path}"
        for entry in sorted(path.iterdir()):
            target = entry.resolve(strict=True)
            assert target.is_relative_to(rtl), f"External library resource: {entry}"
            if target.is_dir():
                copy_directory(entry, ancestors | {resolved})
            else:
                assert target.is_file(), f"Missing regular dependency: {entry}"
                files[entry.relative_to(rtl).as_posix()] = target.read_bytes()

    copy_directory(rtl, set())
    assert files, "Missing delivered RTL"
    return files


@contextmanager
def scratch_library(sources):
    with tempfile.TemporaryDirectory(dir="/tmp") as name:
        parent = Path(name)
        parent.chmod(0o711)
        rtl = parent / "rtl"
        work = parent / "work"
        for directory in (rtl, work):
            directory.mkdir()
            os.chown(directory, 65534, 65534)
        for filename, source in sources.items():
            path = rtl / filename
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(source)
        yield rtl, work


def run_and_reap(command, cwd, env, out, err):
    # This supervisor adopts orphans only from its own tool invocation. Even
    # a child which starts a new session is reparented here when its parent exits.
    supervisor = os.fork()
    if supervisor == 0:
        status = 125
        try:
            if ctypes.CDLL(None).prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
                raise OSError("Cannot track tool descendants")
            try:
                status = subprocess.run(
                    command, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=out,
                    stderr=err, start_new_session=True, **UNPRIVILEGED,
                ).returncode
            finally:
                while True:
                    children = Path("/proc/thread-self/children").read_text().split()
                    if not children:
                        break
                    for child in map(int, children):
                        try:
                            os.kill(child, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        os.waitpid(child, 0)
        except BaseException as exc:
            os.write(err.fileno(), str(exc).encode(errors="replace"))
            status = 125
        os._exit(status if status >= 0 else 128 - status)
    return os.waitstatus_to_exitcode(os.waitpid(supervisor, 0)[1])


def log_excerpt(stream):
    stream.seek(0, os.SEEK_END)
    stream.seek(max(0, stream.tell() - 6000))
    return stream.read().decode(errors="replace")


def completed(stream, marker):
    stream.seek(0)
    target = marker.encode()
    tail = b""
    while chunk := stream.read(65536):
        data = tail + chunk
        if target in data:
            return True
        tail = data[-(len(target) - 1):]
    return False


def run_tool(command, rtl, marker=None):
    env = {**os.environ, "HOME": str(rtl), "TMPDIR": str(rtl.parent / "work"),
           "MPLCONFIGDIR": str(rtl.parent / "work")}
    env.pop("PYTEST_CURRENT_TEST", None)
    with tempfile.TemporaryFile() as out, tempfile.TemporaryFile() as err:
        try:
            status = run_and_reap(command, rtl, env, out, err)
        finally:
            subprocess.run(["find", rtl.parent, "-type", "l", "-delete"], check=True)
        assert status == 0, f"{Path(command[0]).name} exited {status}\n{log_excerpt(out)}\n{log_excerpt(err)}"
        if marker:
            assert completed(out, marker), f"Bench did not complete\n{log_excerpt(out)}\n{log_excerpt(err)}"


def top_files(rtl):
    # The published shell glob skips dotfiles and does not descend into folders.
    names = sorted(glob.glob("*.v", root_dir=rtl), key=os.fsencode)
    assert names, "No compilation units match *.v"
    return ["./" + name for name in names]


@cache
def icarus_headers():
    listing = subprocess.check_output(["/usr/bin/dpkg-query", "-L", "iverilog"], text=True)
    return {path.resolve() for entry in listing.splitlines() if "/ivl/include/" in entry
            and (path := Path(entry)).is_file()}


def validate_dependencies(paths, rtl, permitted, builtin=False):
    for item in paths:
        path = (rtl / item).resolve()
        allowed = path.is_relative_to(rtl) or path.is_relative_to(OUT / "rtl") or path in permitted
        if builtin and path.is_relative_to("/usr/share/yosys"):
            allowed = True
        assert allowed, f"Library requires an external design dependency: {item}"


def icarus(rtl, work, units, bench, top, image=None):
    deps = work / "icarus.dep"
    target = ["-o", image] if image is not None else ["-tnull"]
    run_tool(["/usr/bin/iverilog", "-g2012", "-s", top, f"-Mprefix={deps}",
              *target, *units, bench], rtl)
    assert deps.is_file(), "Icarus did not record its design inputs"
    permitted = {bench, *(Path(p) for p in units if Path(p).is_absolute())}
    paths = []
    for line in deps.read_text().splitlines():
        kind, separator, filename = line.partition(" ")
        assert separator and kind in ("I", "M"), f"Unrecognized dependency entry {line}"
        paths.append(filename)
    validate_dependencies(paths, rtl, permitted | icarus_headers())


def yosys(rtl, work, commands, *permitted):
    deps = work / "yosys.dep"
    run_tool(["/usr/bin/yosys", "-Q", "-E", deps, "-p", commands], rtl)
    assert deps.is_file(), "Yosys did not record its design inputs"
    # Yosys emits Make-escaped file names, including spaces. The targets
    # are the explicitly selected output filenames in the scratch directory.
    _, sep, inputs = deps.read_text().partition(":")
    assert sep, "Malformed synthesis dependency record"
    # This is a dependency record, so apostrophes and quotation marks in
    # filenames are literal characters. Only escaped whitespace joins words.
    words = re.findall(r"(?:\\[ \t]|\S)+", inputs)
    paths = [re.sub(r"\\([ \t])", r"\1", word) for word in words]
    validate_dependencies(paths, rtl, permitted, builtin=True)


@cache
def simulation_monitor():
    directory = tempfile.TemporaryDirectory(dir="/tmp")
    atexit.register(directory.cleanup)
    root = Path(directory.name)
    root.chmod(0o711)
    program = root / "vvp_files"
    subprocess.run(["/usr/bin/cc", "-O2", "-Wall", "-o", program, GRADER / "vvp_files.c"],
                   check=True)
    return program


def simulate_in(rtl, work, units, bench, inspect_image=None):
    case = rtl.parent / "case.sv"
    completion = "BENCH_OK_" + hashlib.sha256(bench.encode()).hexdigest()[:24]
    case.write_text(bench.replace('"BENCH_OK"', f'"{completion}"'))
    image = work / "capture.vvp"
    icarus(rtl, work, units, case, "tb", image)
    assert image.is_file(), "No executable simulation image"
    if inspect_image:
        inspect_image(image)
    run_tool([simulation_monitor(), rtl.parent, OUT / "rtl", "/usr/bin/vvp", image], rtl, marker=completion)


def simulate(sources, bench, inspect_image=None):
    with scratch_library(sources) as (rtl, work):
        simulate_in(rtl, work, top_files(rtl), bench, inspect_image)


@pytest.fixture(params=BUILDS, ids=build_id)
def build(request):
    return request.param


@pytest.fixture(params=DOMAIN_BUILDS, ids=build_id)
def domain_build(request):
    return request.param


def window_bench(build, scenario):
    a, d, g = build
    return "\n".join((
        "`timescale 1ns/1ps", "module tb;",
        f"localparam integer A={a}, D={d}, G={g}, SCENARIO={scenario};",
        (GRADER / "window_core.svh").read_text(),
        (GRADER / "window_stalls.svh").read_text(),
        (GRADER / "window_interactions.svh").read_text(),
        (GRADER / "window_cases.svh").read_text(),
        "endmodule",
    ))


def window(sources, build, scenario):
    simulate(sources, window_bench(build, scenario))


def test_complete_library(sources):
    """Ship the complete library in /app/out/rtl with the original filenames and any new .v dependencies."""
    assert set(LIBRARY) <= sources.keys(), f"Missing library files: {set(LIBRARY) - sources.keys()}"


def legacy_tools():
    spec = importlib.util.spec_from_file_location("legacy_checks", GRADER / "legacy_checks.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_existing_library_callers(sources):
    """Existing integrations can still bind the shipped raw core, FIFO and stream helpers."""
    bench, top, prefix = legacy_tools().protected_callers(INPUTS, sources)
    with scratch_library(sources) as (rtl, work):
        case = rtl.parent / "callers.sv"
        case.write_text(bench)
        icarus(rtl, work, top_files(rtl), case, top)
        metadata = work / "public.json"
        bindings = work / "parameters.txt"
        yosys(rtl, work, f"read_verilog -sv *.v; read_verilog -sv {case}; hierarchy -check -top {top}; "
                         f"tee -q -o {bindings} chparam -list; proc; write_json {metadata}", case, bindings)
        legacy_tools().check_public_signatures(metadata, bindings, prefix)


@pytest.mark.parametrize("module", ["sfifo", "axisrle", "skidbuffer"])
def test_helper_positional_pin_behavior(sources, module):
    """Existing FIFO and stream callers retain positional parameter, port and pin behavior."""
    prefix = legacy_tools().comparison_namespace(sources)
    bench = legacy_tools().baseline_library(INPUTS, prefix) + "\n" + (
        GRADER / "helper_pins.sv").read_text().replace("baseline_", prefix)
    simulate(sources, f"`define COMPARE_{module}\n" + bench)



@pytest.mark.parametrize("client", ["raw", "compressed", "supplied-compressed", "supplied-raw", "long-raw"])
def test_existing_integration_pins(sources, client):
    """“Preserve the existing ports, positional parameter bindings and pin behavior” with the option off."""
    bench, reference = legacy_tools().legacy_bench(INPUTS, sources, client)
    simulate(sources, bench, lambda image: legacy_tools().check_legacy_ports(image, reference))


@pytest.mark.parametrize("client", ["synthesized-raw", "synthesized-compressed", "synthesized-positional"])
def test_synthesized_default_pins(sources, client):
    """The synthesized raw and compressed clients preserve the supplied default-mode pin behavior."""
    bench, _ = legacy_tools().legacy_bench(INPUTS, sources, client)
    module = "memscopc" if client == "synthesized-compressed" else "memscope"
    with scratch_library(sources) as (rtl, work):
        netlist = work / "hardware.v"
        top = module
        configure = ""
        allowed = []
        if client == "synthesized-positional":
            top = legacy_tools().comparison_namespace(sources) + "bound"
            caller = rtl.parent / "binding.sv"
            caller.write_text(legacy_tools().positional_raw_caller(INPUTS, top))
            allowed.append(caller)
            configure = (f"read_verilog -sv {caller}; chparam -set C_AXI_ADDR_WIDTH 16 "
                         "-set C_AXI_DATA_WIDTH 64 -set C_AXI_ID_WIDTH 3 "
                         "-set OPT_TREADY_WHILE_IDLE 0 -set LGFIFO 5 -set LGLEN 15 "
                         f"-set AXI_ID 5 -set HOLDOFFBITS 12 -set DEF_HOLDOFF 37 {top}; ")
        rename = f"rename -top {module}; " if top != module else ""
        yosys(rtl, work, f"read_verilog -sv *.v; {configure}synth -top {top} -flatten -run begin:fine; "
                         f"techmap; opt -fast; check -assert; {rename}write_verilog -noattr {netlist}",
              *allowed)
        simulate_in(rtl, work, [netlist], bench)


def test_staging_registers_and_decoupled_lite(sources, build):
    """“Writes merge by WSTRB and reads return the staged values” even with skewed AW and W requests."""
    window(sources, build, 1)


def test_unsigned_descriptor_domain(sources, domain_build):
    """“B is w-aligned, N>0, B+N*w<=2^A in mathematical integers and H<N.”"""
    window(sources, domain_build, 2)


def test_continuous_one_slot_ring(sources, domain_build):
    """“Accept one sample every clock from the first acquisition edge through closure” also applies to N=1."""
    window(sources, domain_build, 3)


def test_burst_sizing_and_full_rate(sources, domain_build):
    """“At byte address p use a full-width, all-lanes INCR burst” subject to the FIFO, ring and page limits."""
    window(sources, domain_build, 4)


def test_trigger_qualification_and_holdoff(sources, build):
    """“Include exactly H subsequent accepted samples” after a qualifying beat sees its pre-edge flags."""
    window(sources, build, 5)


def test_reassign_unprimed_prefix(sources, domain_build):
    """“On that edge a stream handshake still belongs to the outgoing capture and sees the prior controls.”"""
    window(sources, domain_build, 6)


def test_reassign_triggered_capture(sources, build):
    """ARM “drains its entire accepted prefix to its own ring” across independently stalled AW, W and B."""
    window(sources, build, 7)


def test_reassign_during_drain_suppresses_interrupt(sources, build):
    """“A valid ARM closes that capture, suppresses its pending completion interrupt” even during drain."""
    window(sources, build, 8)


def test_stopped_chronology_lanes_and_rewind(sources, domain_build):
    """“Return the last N samples oldest first, low 32-bit lane first and cyclically” including after rewind."""
    window(sources, domain_build, 9)


def test_reassign_with_old_read(sources, build):
    """A Lite read can remain outstanding across ARM while outgoing DDR reads must retire."""
    window(sources, build, 10)


def test_lite_read_preupdate_snapshot(sources, build):
    """“Read registers as of the AR handshake, before same-edge updates.”"""
    window(sources, build, 11)


def test_invalid_arm_has_no_effect_on_acquisition(sources, build):
    """“An invalid ARM returns SLVERR and has no effect” during an active and a triggered capture."""
    window(sources, build, 12)


def test_capture_with_independent_backpressure(sources, build):
    """“Write each accepted sample i once, in stream order” when the AXI4 channels stall independently."""
    window(sources, build, 13)


def test_valid_with_waiting_destinations(sources, build):
    """AXI sources offer responses, write addresses, write data and read addresses to waiting peers."""
    window(sources, build, 14)


def test_live_capture_counter_widths(sources, domain_build):
    """The admitted sample capacity and holdoff govern burst limits, priming and the final accepted beat."""
    window(sources, domain_build, 15)


def test_trigger_edges_and_status_snapshots(sources, build):
    """Trigger only on an eligible accepted beat and snapshot status before each same-edge transition."""
    window(sources, build, 16)


@pytest.mark.parametrize("scenario", [17, 18], ids=["stalled-address", "withheld-data"])
def test_stopped_read_overlaps_arm_execution(sources, build, scenario):
    """A stopped-data request can overlap ARM, with the old allocation fenced through DDR retirement."""
    window(sources, build, scenario)


def test_arm_command_while_ddr_read_response_withheld(sources, build):
    """Accept the Lite command independently of DDR progress while fencing its acknowledgment on R."""
    window(sources, build, 19)


def test_arm_control_ignores_unassigned_fields(sources, build):
    """ARM accepts exactly the descriptor predicate, including status bits written back by firmware."""
    window(sources, build, 20)


def test_capacity_and_holdoff_at_counter_boundaries(sources, domain_build):
    """Prime after N samples and include H subsequent transfers using the full active descriptor."""
    window(sources, domain_build, 21)


def test_trigger_on_backpressured_word(sources, build):
    """An external trigger qualifies only on the accepted sampling edge."""
    window(sources, build, 22)


def test_one_slot_response_retirement(sources, build):
    """A closed capture fences all outstanding B handshakes before interrupt or ARM acknowledgment."""
    window(sources, build, 23)


def test_arm_on_natural_closing_edge(sources, build):
    """ARM retains its coincident outgoing sample and suppresses the capture's pending completion."""
    window(sources, build, 24)


def test_rejected_reassignment_preserves_pending_completion(sources, build):
    """Invalid ARM has no effect on the pending interrupt, retirement or stopped image."""
    window(sources, build, 25)


def test_register_command_acceptance_during_ddr_stalls(sources, build):
    """Accept Lite AW and W independently of DDR progress for the in-scope register writes."""
    window(sources, build, 26)


def test_open_burst_length_during_stream_gaps(sources, domain_build):
    """Only the last burst of a closed capture may shorten below the prescribed length."""
    window(sources, domain_build, 27)


def synthesize_window(sources, build):
    a, d, g = build
    configure = (
        "read_verilog -sv *.v; chparam -set OPT_WINDOWED 1 -set OPT_TREADY_WHILE_IDLE 0 "
        f"-set C_AXI_ADDR_WIDTH {a} -set C_AXI_DATA_WIDTH {d} -set LGFIFO {g} memscope; "
        "hierarchy -check -top memscope"
    )
    full_synthesis = "synth -top memscope; check -assert"
    with scratch_library(sources) as (rtl, work):
        netlist = work / "hardware.sv"
        # The standard synthesis checks and the simulated RAM-preserving
        # circuit start from the same configured hardware. Keep word-level
        # arithmetic for simulation and inspect bit-level logic for AXI paths.
        connectivity = work / "connectivity.json"
        commands = (f"{configure}; design -push-copy; {full_synthesis}; design -pop; "
                    "synth -top memscope -flatten -noalumacc -run begin:fine; opt -full; "
                    f"check -assert; write_verilog -noattr {netlist}; "
                    "setattr -mod -unset keep_hierarchy; setattr -unset keep_hierarchy; "
                    "flatten -wb -noscopeinfo; techmap; opt -fast; "
                    f"write_json {connectivity}")
        yosys(rtl, work, commands)
        assert netlist.is_file() and connectivity.is_file(), "Missing synthesized circuit"
        return netlist.read_bytes(), json.loads(connectivity.read_bytes())["modules"]["memscope"]


@pytest.fixture(scope="module")
def hardware(sources):
    @cache
    def for_build(build):
        return synthesize_window(sources, build)
    return for_build


def synthesized_window(hardware, sources, build, scenario):
    circuit, _ = hardware(build)
    with scratch_library(sources) as (rtl, work):
        netlist = rtl.parent / "hardware.v"
        netlist.write_bytes(circuit)
        simulate_in(rtl, work, [netlist], "`define FIXED_WINDOW_NETLIST\n" + window_bench(build, scenario))


@pytest.mark.parametrize(("build", "scenario"), [
    pytest.param(build, scenario, id=f"{build_id(build)}-{scenario}")
    for build in DOMAIN_BUILDS
    for scenario in (range(1, 28) if build in BUILDS else WIDTH_SCENARIOS)
])
def test_synthesized_window_transactions(hardware, sources, build, scenario):
    """The synthesized window implements the same stream, DDR and firmware contract at its pins."""
    synthesized_window(hardware, sources, build, scenario)


def test_no_combinational_axi_interface_paths(hardware, domain_build):
    """AXI interfaces have no combinational paths from their channel inputs to their outputs."""
    _, module = hardware(domain_build)
    loads = defaultdict(set)
    for cell in module["cells"].values():
        kind = cell["type"]
        pins = cell["connections"]
        if kind == "$mem_v2":
            # Each asynchronous read port propagates only its own address.
            p = {k: int(v, 2) for k, v in cell["parameters"].items()
                 if k in ("RD_CLK_ENABLE", "RD_PORTS", "ABITS", "WIDTH")}
            for port in range(p["RD_PORTS"]):
                if (p["RD_CLK_ENABLE"] >> port) & 1:
                    continue
                addr_start = port * p["ABITS"]
                data_start = port * p["WIDTH"]
                data_bits = pins["RD_DATA"][data_start:data_start + p["WIDTH"]]
                for bit in pins["RD_ADDR"][addr_start:addr_start + p["ABITS"]]:
                    loads[bit].update(b for b in data_bits if isinstance(b, int))
            continue
        if kind.startswith("$_") and "DFF" in kind:
            if not kind.startswith("$_SDFF"):
                for control in ("R", "S", "AD", "L"):
                    for bit in pins.get(control, []):
                        loads[bit].update(b for b in pins["Q"] if isinstance(b, int))
            continue
        # Transparent latches can propagate both data and gate changes
        # between clock edges, as can asynchronous set/reset latches.
        if not kind.startswith("$_"):
            continue
        outputs = {b for name, bits in pins.items() if cell["port_directions"][name] == "output"
                   for b in bits if isinstance(b, int)}
        for name, bits in pins.items():
            if cell["port_directions"][name] == "input":
                for bit in bits:
                    loads[bit].update(outputs)
    for prefix in ("S_AXIL_", "M_AXI_", "S_AXIS_"):
        ports = {n: v for n, v in module["ports"].items() if n.startswith(prefix)}
        reached = {b for v in ports.values() if v["direction"] == "input"
                   for b in v["bits"] if isinstance(b, int)}
        todo = list(reached)
        while todo:
            new = loads[todo.pop()] - reached
            reached.update(new)
            todo.extend(new)
        for name, port in ports.items():
            if port["direction"] == "output":
                assert reached.isdisjoint(port["bits"]), f"Combinational path from {prefix} inputs to {name}"


@pytest.mark.parametrize("build", (*BUILDS, (12, 32, 2)), ids=build_id)
@pytest.mark.parametrize("circuit", ("rtl", "synthesized"))
def test_stage_next_allocation_while_arm_is_outstanding(sources, hardware, build, circuit):
    """A staging write can overlap ARM while the accepted descriptor and response order are preserved."""
    if circuit == "rtl":
        window(sources, build, 28)
    else:
        synthesized_window(hardware, sources, build, 28)


@pytest.mark.parametrize("build", (*BUILDS, (12, 32, 2)), ids=build_id)
@pytest.mark.parametrize("circuit", ("rtl", "synthesized"))
def test_stage_both_registers_behind_arm(sources, hardware, build, circuit):
    """While ARM awaits B, both staging registers accept writes in either order for the following ARM."""
    if circuit == "rtl":
        window(sources, build, 29)
    else:
        synthesized_window(hardware, sources, build, 29)


@pytest.mark.parametrize("build", (*BUILDS, (12, 32, 2)), ids=build_id)
@pytest.mark.parametrize("circuit", ("rtl", "synthesized"))
def test_two_outstanding_reads_across_capture_handoff(sources, hardware, build, circuit):
    """Accept the second AR under stalls, return snapshots in AR order, and fence both old fetches at ARM."""
    if circuit == "rtl":
        window(sources, build, 30)
    else:
        synthesized_window(hardware, sources, build, 30)
