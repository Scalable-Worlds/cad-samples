"""Regenerate the native comparator and read the values displayed in its stock tables."""

from collections import Counter, defaultdict
from functools import lru_cache
from dataclasses import dataclass
import base64
import binascii
from pathlib import Path
from urllib.parse import unquote, unquote_to_bytes, urlsplit
from xml.etree import ElementTree
from decimal import Decimal
import math
import re
import unicodedata


def build(preset, angle, width, cut):
    """Use the supplied package on the trusted import path for one native build."""
    import FreeCAD
    from openafpm_cad_core.dimension_tables import get_dimension_tables
    from openafpm_cad_core.get_default_parameters import get_default_parameters
    from openafpm_cad_core.high_end_stop_cells import high_end_stop_cells
    from openafpm_cad_core.load import load_turbine

    assert math.isfinite(angle) and math.isfinite(width), "Stop settings must be finite"
    config = get_default_parameters(preset)
    nominal_boom = config["furling"]["BoomLength"]
    config = {**config, "furling": dict(config["furling"])}
    if cut is not None:
        config["furling"]["BoomLength"] = cut
    inputs = [config[group] for group in ("magnafpm", "furling", "user")]
    # The supplied spreadsheet generator uses degree quantities for the rotation
    # and plain millimetre magnitudes for lengths. Replace just the two formulas.
    values = {"MaximumFurlAngle": f"{angle:.17g} deg", "HighEndStopWidth": f"{width:.17g}"}
    edits = [(cell, cell.content) for row in high_end_stop_cells for cell in row
             if cell.alias in values]
    try:
        for cell, original in edits:
            cell.content = values[cell.alias]
        wind, master = load_turbine(*inputs)
    finally:
        for cell, original in edits:
            cell.content = original
    tables = get_dimension_tables(master, FreeCAD.getDocument("Alternator"))
    return wind, master, tables, config, nominal_boom


def _nodes(element, tag):
    if element.get("tagName", "").casefold() == tag:
        yield element
    else:
        for child in element.get("children", []):
            yield from _nodes(child, tag)


def _text(element):
    value = element.get("properties", {}).get("textContent")
    if value is not None:
        return str(value).strip()
    return " ".join(filter(None, map(_text, element.get("children", [])))).strip()


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
SCALES = {
    "μm": ("length", Decimal(".001")), "um": ("length", Decimal(".001")),
    "mm": ("length", Decimal(1)), "cm": ("length", Decimal(10)),
    "m": ("length", Decimal(1000)), "in": ("length", Decimal("25.4")),
    "kg": ("mass", Decimal(1)), "g": ("mass", Decimal(".001")),
    "rad": ("angle", Decimal(str(math.degrees(1)))), "deg": ("angle", Decimal(1)), "°": ("angle", Decimal(1)),
}
UNIT = r"(?:[μu]m\b|mm\b|cm\b|kg\b|m\b|g\b|in\b|rad\b|deg\b|°)"
QUANTITY = re.compile(fr"({NUMBER})\s*({UNIT})?")
STOCK_PAIR = re.compile(fr"(?<![a-z])(?P<thread>m\s*)?(?P<dia>{NUMBER})\s*"
                        fr"(?P<du>{UNIT})?\s*[×x]\s*(?P<length>{NUMBER})\s*(?P<lu>{UNIT})?")
COUNT = r"(?:pieces?|pcs?\.?|each|ea\.?)\b"


@lru_cache(maxsize=4096)
def _clean(value):
    # Normalize superscript scientific notation before compatibility folding.
    text = re.sub(r"10([⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻]+)", lambda m: "10^" + unicodedata.normalize("NFKC", m[1]), str(value))
    value = unicodedata.normalize("NFKC", text).lower().replace("−", "-")
    value = re.sub(fr"({NUMBER})\s*[×x·]\s*10\s*(?:\^|\*\*)\s*([+-]?\d+)", r"\1e\2", value)
    for pattern, symbol in (
        (r"\bmicromet(?:er|re)s?\b|\bmicrons?\b", "μm"),
        (r"\bmillimet(?:er|re)s?\b", "mm"), (r"\bcentimet(?:er|re)s?\b", "cm"),
        (r"\bmet(?:er|re)s?\b", "m"), (r"\bkilograms?\b", "kg"),
        (r"\bradians?\b", "rad"), (r"\bgrams?\b", "g"), (r"\bdegrees?\b", "deg"), (r"\binch(?:es)?\b", "in"),
        (r"\b(?:outside|external)\s+dia(?:meter)?\b|\bo\.?d\.?\b", "outer diameter"),
        (r"\binside\s+dia(?:meter)?\b|\bi\.?d\.?\b", "inner diameter"),
        (r"\bdia\b|[ø⌀]", " diameter "), (r"\bmold\b", "mould"),
        (r"\bcentre\b", "center"), (r"\bdisc\b", "disk"),
        (r"\byaw(?: bearing)? (?:pipe|tube)\b", "yaw pipe"), (r"\bpitch circle diameter\s*\(pcd\)|\bpcd\b", "pitch circle diameter"),
        (r"\bmass\b", "weight"), (r"\bqty\b|\bquantity\b", "number"),
        (r"\b(?:section|figure|fig\.?|page)\s*#?\s*\d+\b", ""),
    ):
        value = re.sub(pattern, symbol, value)
    return value


@lru_cache(maxsize=4096)
def _name(value):
    """Normalize engineering terms and retain dimensioned stock identifiers."""
    def dimension(match):
        number = Decimal(match[1])
        kind, factor = SCALES.get(match[2], ("number", Decimal(1)))
        return " " + kind + format((number * factor).normalize(), "f") + " "
    value = re.sub(fr"(?:\(|\[)\s*{UNIT}\s*(?:\)|\])\s*$", "", _clean(value))
    # Parenthesized drafting symbols qualify a written-out quantity.
    value = re.sub(r"(?<=\w)\s*\([a-z]\)\s*$", "", value)
    value = QUANTITY.sub(dimension, value)
    ignored = {"for", "the", "at", "of", "to", "from", "and", "by", "per", "in",
               "nominal", "required", "sizes", "size", "dimensions", "dimension",
               "estimated", "mm", "cm", "kg", "deg"}
    return tuple(sorted(w[:-1] if w in {"bolts", "nuts", "studs", "screws", "washers",
                                        "magnets", "coils", "pipes", "disks", "pieces"} else w
                        for w in re.findall(r"[a-z]+(?:[-+]?\d+(?:\.\d+)?)?|\d+", value)
                        if w not in ignored))


def _unit_context(label):
    match = re.search(fr"(?:\(|\[)\s*({UNIT})\s*(?:\)|\])\s*$", _clean(label))
    return match[1] if match else None


@dataclass
class Table:
    title: str
    rows: list
    notes: str


def _svg_element(node):
    props = dict(node.get("properties", {}))
    text = str(props.pop("textContent", ""))
    xml = ElementTree.Element(node["tagName"], {k: str(v) for k, v in props.items()})
    xml.text = text
    for child in node.get("children", []):
        xml.append(_svg_element(child))
    return xml


def _image_sources(element):
    for node in _nodes(element, "img"):
        yield node.get("properties", {}).get("src", "")
    for tag, prop in (("object", "data"), ("embed", "src")):
        for node in _nodes(element, tag):
            attrs = node.get("properties", {})
            if attrs.get("type", "image/").startswith("image/"):
                yield attrs.get(prop, "")
    for node in _nodes(element, "svg"):
        yield ElementTree.tostring(_svg_element(node), encoding="utf-8")


def _tables(elements):
    for table in (node for element in elements for node in _nodes(element, "table")):
        rows = list(_nodes(table, "tr"))
        if not rows:
            continue
        headings = list(_nodes(table, "caption")) or list(_nodes(table, "thead"))
        title = _text(headings[0]) if headings else _text(rows[0])
        groups = list(_nodes(table, "tbody"))
        body = [r for g in groups for r in _nodes(g, "tr")] if groups else rows
        cells = [[_text(c) for c in r.get("children", [])
                  if c.get("tagName", "").lower() in {"td", "th"}] for r in body]
        headers = [[_text(c) for c in r.get("children", [])] for g in _nodes(table, "thead")
                   for r in _nodes(g, "tr")]
        wide = next((h for h in reversed(headers) if len(h) >= 2), None)
        if not wide and len(cells) >= 2 and len(cells[0]) > 1 and all(
                not QUANTITY.search(c) for c in cells[0]) and all(
                len(c) == len(cells[0]) and all(QUANTITY.search(x) for x in c) for c in cells[1:]):
            wide, cells = cells[0], cells[1:]
        facts = []
        for row in cells:
            if wide and len(wide) == len(row):
                facts.extend((label, value) for label, value in zip(wide, row) if label and value)
            elif len(row) >= 2 and row[0] and any(row[1:]):
                facts.append((row[0], " ".join(c for c in row[1:] if c)))
        yield Table(title, facts, " ".join(_text(g) for g in _nodes(table, "tfoot")))


def _scopes(native, returned):
    """Identify scope from headings or the subjects of rows, independent of values."""
    headings = [_name(t.title) for t in native]
    keys = [set(_name(label) for label, _ in t.rows) for t in native]
    result = []
    for table in returned:
        heading = _name(table.title)
        ids = {_name(label) for label, _ in table.rows}
        scores = [(heading == h, len(ids & names),
                   len(set(heading) & set(h))) for h, names in zip(headings, keys)]
        top = max(scores)
        result.append([i for i, score in enumerate(scores) if score == top])
    return result


def _reading(value, inherited=None, hardware=False):
    """Retain diameter, thread, length and count roles when cells reorder quantities."""
    value = _clean(value).replace("—", " ").replace("–", " ")
    out = defaultdict(list)

    def append(raw, role, unit=None):
        n = Decimal(raw)
        _, factor = SCALES.get(unit, (role, Decimal(1)))
        # Thread and bare hardware diameters identify stock. Counts are exact.
        discrete = role in {"thread", "count", "number", "diameter"}
        precision = Decimal(0) if discrete else Decimal(".5").scaleb(n.as_tuple().exponent) * factor
        out[role].append((n * factor, precision))

    def pair(match):
        dia_unit = match['du'] or (None if match['thread'] else match['lu'])
        append(match['dia'], "thread" if match['thread'] else "diameter", dia_unit)
        append(match['length'], "length", match['lu'] or match['du'] or "mm")
        return " "
    if hardware:
        # In "4 × M12" and "4 × 5 × 30 mm" the leading multiplier is
        # the piece count. The remaining product designates the stock.
        prefix = fr"({NUMBER})\s*[x×]\s*(?=(?:m\s*{NUMBER}|{NUMBER}\s*(?:{UNIT})?\s*[x×]))"
        def quantity(match):
            append(match[1], "count")
            return " "
        value = re.sub(prefix, quantity, value)
    value = STOCK_PAIR.sub(pair, value)

    def apply(pattern, role, field):
        def record(match):
            append(match[field], role)
            return " "
        return re.sub(pattern, record, value)

    value = apply(fr"\bm\s*({NUMBER})", "thread", 1)
    value = apply(fr"({NUMBER})\s*{COUNT}", "count", 1)
    value = apply(fr"(?:\bnumber\b|\bpcs\b|\bpieces\b)\s*[:=]?\s*({NUMBER})", "count", 1)
    if hardware:
        value = apply(fr"({NUMBER})\s*(?:bolts?|screws?|nuts?|washers?|studs?|rods?)\b", "count", 1)
        # Some tables label the two stock dimensions instead of using A × B.
        for pattern in (fr"diameter\s*[:=]?\s*({NUMBER})\s*({UNIT})?",
                        fr"({NUMBER})\s*({UNIT})?\s*diameter"):
            def diameter(match):
                append(match[1], "diameter", match[2] or inherited)
                return " "
            value = re.sub(pattern, diameter, value)
    for match in QUANTITY.finditer(value):
        raw, unit = match.groups()
        unit = unit or inherited
        kind = SCALES[unit][0] if unit else "number"
        append(raw, kind, unit)
    return dict(out)


def _agrees(label, expected, actual_label, actual, title):
    left = _reading(expected, _unit_context(label))
    right = _reading(actual, _unit_context(actual_label) or _unit_context(title),
                     hardware=bool({"diameter", "thread"} & left.keys()))
    # Material grades include qualitative designations as well as numbers.
    if not left or "material" in _name(label):
        return _name(expected) == _name(actual)
    if set(left) == {"number"} and set(right) == {"count"}:
        right = {"number": right["count"]}
    if left.keys() != right.keys():
        return False
    for role, values in left.items():
        others = right[role]
        if not all(any(abs(a - b) <= tol for a, tol in values) for b, _ in others):
            return False
        if not all(any(abs(a - b) <= tol for b, _ in others) for a, tol in values):
            return False
    return True


def table_differences(expected_tables, returned_tables):
    """Match stock and quantities in context, preserving native row multiplicity."""
    native, actual = list(_tables(expected_tables)), list(_tables(returned_tables))
    keys = defaultdict(list)
    for index, table in enumerate(native):
        for label, value in table.rows:
            keys[_name(label)].append((index, label, value))
    seen = Counter()
    errors = []
    owners_by_table = _scopes(native, actual)
    for table, scopes in zip(actual, owners_by_table):
        for label, value in table.rows:
            label_terms = _name(label)
            candidates = keys.get(label_terms, [])
            if not candidates:
                # A component may be identified by either its row or its table.
                # Compare the residual subject only within that table context.
                for i in scopes:
                    context = set(_name(native[i].title))
                    subject = set(label_terms) - context - {"finished", "stock", "bearing"}
                    if not subject:
                        continue
                    candidates = candidates + [(i, key, val) for key, val in native[i].rows
                                                if set(_name(key)) - context - {"finished", "stock", "bearing"} == subject]
            scoped = [row for row in candidates if row[0] in scopes]
            # A known component heading assigns a manufacturing context to
            # all its rows. Consolidated schedules may supply their own scope.
            named_table = any(_name(table.title) == _name(native[i].title) for i in scopes)
            options = scoped if named_table else scoped or candidates
            if candidates and not options:
                errors.append(f"{table.title} / {label} is assigned to the wrong component")
                continue
            if not options:
                # An explanatory row is allowed. A newly asserted stock size in
                # the native stock-total schedule is a new manufacturing claim.
                stock_subject = [w for w in _name(label) if re.match(r"length\d", w)]
                if stock_subject and "diameter" in _name(label):
                    errors.append(f"Unexpected stock specification {label}: {value}")
                continue
            matches = [key for key in options
                       if _agrees(key[1], key[2], label, value, table.title)]
            if matches:
                key = min(matches, key=lambda k: seen[k])
                seen[key] += 1
            else:
                errors.append(f"{table.title} / {label}: {value} disagrees with the native specification")
    # A changed alloy is a substantive change to the shop instructions. Match
    # the material in component context, independent of where the note is laid out.
    for table in native:
        if "stainless steel" in _clean(table.notes):
            scope = native.index(table)
            local = " ".join(" ".join([t.notes, *(v for _, v in t.rows)])
                             for t, owners in zip(actual, owners_by_table) if scope in owners)
            if re.search(r"\b(?:are|material\s*[:=]|made of)\s+(?:uncoated\s+)?(?:mild|carbon) steel\b",
                         _clean(local)):
                errors.append(f"Conflicting fastener material in {table.title}")
    for rows in keys.values():
        for key, count in Counter(rows).items():
            if seen[key] < count:
                errors.append(f"{native[key[0]].title} / {key[1]}: expected {count} matching row(s), got {seen[key]}")
    return errors


def _resource(src, package, bases, visited):
    """Resolve an illustration and local resources referenced by an SVG."""
    if isinstance(src, bytes):
        data = src
    else:
        if not isinstance(src, str) or not src:
            return False
        path = urlsplit(src)
        if not path.path and path.fragment:
            return True
        if path.scheme == "data":
            header, _, payload = src.partition(",")
            if not header.lower().startswith("data:image/"):
                return False
            data = unquote_to_bytes(payload)
            if ";base64" in header.lower():
                data = base64.b64decode(data, validate=True)
        elif path.scheme in ("", "file") and not path.netloc:
            for base in bases:
                candidate = (base / unquote(path.path)).resolve()
                if candidate.is_relative_to(package) and candidate.is_file():
                    if candidate in visited:
                        return True
                    visited.add(candidate)
                    data = candidate.read_bytes()
                    bases = [candidate.parent]
                    break
            else:
                return False
        else:
            return False
    if data.removeprefix(b"\xef\xbb\xbf").lstrip().startswith(b"<"):
        svg = ElementTree.fromstring(data)
        return all(_resource(value, package, bases, visited) for elem in svg.iter()
                   if elem.tag.rsplit("}", 1)[-1] in {"image", "use"}
                   for attr, value in elem.attrib.items() if attr.rsplit("}", 1)[-1] == "href")
    return bool(data)


def illustration_inventory(tables, package):
    package = Path(package).resolve()
    missing = []
    for element in tables:
        for src in _image_sources(element):
            try:
                valid = _resource(src, package, [package, package / "img", package.parent], set())
            except (OSError, ValueError, SyntaxError, ElementTree.ParseError, binascii.Error):
                valid = False
            if not valid:
                missing.append(str(src)[:160])
    return missing
