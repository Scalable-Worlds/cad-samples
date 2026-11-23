import hashlib
import json
import re
from pathlib import Path


def check_legacy_ports(image, reference_module):
    """Compare Icarus port metadata for each elaborated pair in the legacy bench."""
    # The instance denotes a scope, whose metadata names the external ports
    # even for an explicit mapping such as .S_AXI_ACLK(capture_clock).
    scope = re.compile(
        r'^\S+ \.scope module, "(dut|ref_dut)" "([^"]+)" .*?, (\S+);$'
    )
    port = re.compile(r'^\s+\.port_info (\d+) /(INPUT|OUTPUT|INOUT) (\d+) "([^"]+)";$')
    scopes = {}
    current = None
    with image.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if " .scope " in line:
                current = None
                match = scope.match(line)
                if match:
                    name, module, parent = match.groups()
                    current = []
                    scopes.setdefault(parent, {})[name] = (module, current)
            elif current is not None:
                match = port.match(line)
                if match:
                    current.append(match.groups())
    pairs = [
        instances for instances in scopes.values()
        if "ref_dut" in instances and instances["ref_dut"][0] == reference_module
    ]
    assert pairs, "Missing elaborated original client interface"
    for pair in pairs:
        expected = pair["ref_dut"][1]
        assert expected and "dut" in pair, "Missing elaborated legacy client or its interface"
        actual = pair["dut"][1]
        assert actual == expected, f"Changed public interface of {pair['dut'][0]}: {actual} vs {expected}"


PUBLIC_MODULES = ("memscope", "memscopc", "axisrle", "skidbuffer", "sfifo")


def comparison_namespace(sources):
    identity = hashlib.sha256(b"".join(sources.values())).hexdigest()[:20]
    return f"comparison_{identity}_"


def baseline_library(inputs, namespace):
    names = re.compile(r"\b(" + "|".join(PUBLIC_MODULES) + r")\b")
    return "\n".join(
        names.sub(lambda match: namespace + match[0], (inputs / f"{name}.v").read_text())
        for name in PUBLIC_MODULES
    )


def rename_declaration(inputs, original, new_name):
    """Keep every instantiation in the protected caller bound to the delivery."""
    source = (inputs / f"{original}.v").read_text()
    declaration = rf"(\bmodule\s+){original}(?=\s*#\s*\()"
    renamed, count = re.subn(declaration, lambda m: m[1] + new_name, source)
    assert count == 1, f"Missing declaration of protected caller {original}"
    return renamed


def protected_callers(inputs, sources):
    """Elaborate independent existing callers and every public library interface."""
    prefix = comparison_namespace(sources)
    modules = [baseline_library(inputs, prefix)]
    names = list(PUBLIC_MODULES) + [prefix + name for name in PUBLIC_MODULES]
    for original in ("memscopc", "memscope"):
        name = prefix + "caller_" + original
        modules.append(rename_declaration(inputs, original, name))
        names.append(name)
    top = prefix + "top"
    instances = [f"  {name} caller{i}();" for i, name in enumerate(names)]
    wrapper = "\n".join([f"module {top};", *instances, "endmodule", ""])
    return "\n".join([*modules, wrapper]), top, prefix


def check_public_signatures(metadata, bindings, prefix):
    """Compare external port types and verify independent named parameter bindings."""
    assert metadata.is_file(), "Missing elaborated public interfaces"
    modules = json.loads(metadata.read_bytes())["modules"]
    assert bindings.is_file(), "Missing positional parameter metadata"
    ordered = {}
    current = None
    for line in bindings.read_text().splitlines():
        if line.endswith(":"):
            current = ordered.setdefault(line[:-1], [])
        elif current is not None and line.startswith("  "):
            current.append(line.strip())

    def signature(info):
        return info["direction"], len(info["bits"]), bool(info.get("signed"))

    for name in PUBLIC_MODULES:
        assert name in modules, f"Missing public module {name}"
        shipped = modules[name]
        original = modules[prefix + name]
        assert list(shipped["ports"]) == list(original["ports"]), f"Public port order changed on {name}"
        for port, expected in original["ports"].items():
            assert port in shipped["ports"], f"Missing {name}.{port}"
            actual = shipped["ports"][port]
            assert signature(actual) == signature(expected), f"Public port type changed at {name}.{port}"
        expected_order = ordered.get(prefix + name, [])
        assert expected_order and ordered.get(name, [])[:len(expected_order)] == expected_order, (
            f"Positional parameter bindings changed on {name}: expected {expected_order}")
        for parameter in expected_order:
            actual = shipped["parameter_default_values"][parameter]
            expected = original["parameter_default_values"][parameter]
            assert int(actual, 2) == int(expected, 2), f"Default changed on {name}.{parameter}"


def legacy_bench(inputs, sources, client):
    prefix = comparison_namespace(sources)
    compressed = client in ("compressed", "supplied-compressed", "synthesized-compressed")
    parts = [baseline_library(inputs, prefix), "`undef COMPRESSED_CLIENT"]
    if compressed:
        parts.append("`define COMPRESSED_CLIENT")
    if client.startswith("synthesized-"):
        parts.append("`define FIXED_LEGACY_NETLIST")
    if client == "synthesized-positional":
        parts.append("`define POSITIONAL_LEGACY_NETLIST")
    if client == "long-raw":
        parts.append("`define LONG_LEGACY_CAPTURE")
    driver = Path(__file__).with_name("legacy.sv").read_text().replace("baseline_", prefix)
    driver = driver.replace("compare_legacy", prefix + "client")
    module = "memscopc" if compressed else "memscope"
    if client in ("supplied-compressed", "supplied-raw"):
        independent = prefix + "caller_" + module
        parts.append(rename_declaration(inputs, module, independent))
        driver = driver.replace(f"`define SUBJECT {module}\n", f"`define SUBJECT {independent}\n")
        if client == "supplied-raw":
            parts.append("`define SUPPLIED_RAW_CLIENT")
    reference = prefix + module
    return "\n".join([*parts, driver]), reference

def positional_raw_caller(inputs, top):
    """Use the protected interface to bind the delivered raw circuit by parameter position."""
    source = (inputs / "memscope.v").read_text()
    declaration = re.search(r"^module\s+memscope\s*#\s*\(.*?^\s*\);", source, re.M | re.S)
    assert declaration, "Missing supplied raw interface"
    header = declaration[0].replace("memscope", top, 1)
    return header + """
    memscope #(C_AXI_ADDR_WIDTH, C_AXI_DATA_WIDTH, C_AXI_ID_WIDTH, OPT_TREADY_WHILE_IDLE,
               LGFIFO, LGLEN, AXI_ID, HOLDOFFBITS, DEF_HOLDOFF) bound(.*);
endmodule
"""
