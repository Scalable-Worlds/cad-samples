"""Read native preset and shop quantities with their dimensions."""
from collections import defaultdict
import math
import re


def _key(text):
    return "".join(ch.lower() for ch in str(text) if ch.isalnum())


def _scalar(value):
    if hasattr(value, "Unit") and hasattr(value, "Value"):
        return {"unit": value.Unit.Type, "value": value.Value}
    return value if isinstance(value, (str, int, float)) else None


def _published_rows(obj, cells):
    """Extract explicit component cut rows, and the optional clearance input."""
    rows = defaultdict(dict)
    for cell, value in cells.items():
        col, row = re.fullmatch(r"([A-Z]+)([1-9]\d*)", cell).groups()
        rows[int(row)][col] = value
    result = []
    def shop_context(text):
        return bool(re.search(r"cut\s*(?:list|length)|schedule|fabrication|stock", text, re.I))

    schedule = shop_context(obj.Label)
    length_columns = set()
    for number, row in sorted(rows.items()):
        if len(row) == 1:
            heading = next(iter(row.values()))
            if isinstance(heading, str) and shop_context(heading):
                schedule = True
        headers = {col for col, label in row.items() if isinstance(label, str) and
                   _key(label) in {"cutlength", "cutlengthmm", "lengthmm", "finishedlength", "length"}}
        if headers:
            length_columns = headers
        for col, label in row.items():
            if not isinstance(label, str):
                continue
            caption = label
            label = _key(caption)
            kind = None
            if label == "vaneclearancefraction":
                kind = "fraction"
            elif "yaw" in label and ("pipe" in label or "tube" in label):
                # Totals by diameter include other components and are a
                # different shop quantity from the cut length of this pipe.
                if any(word in label for word in ("diameter", "radius", "total", "thickness", "projected")):
                    continue
                # A model calculation can have a local YawPipeLength input
                # (e.g. a projected stock length). A cut claim identifies a
                # shop operation, either explicitly or through its schedule.
                if ("length" in label and any(word in label for word in ("cut", "finished", "stock"))) or (
                        (schedule or length_columns) and ("length" in label or
                        label in {"yawpipe", "yawbearingpipe", "yawbearingtube", "yawtube"})):
                    kind = "yaw_cut"
            if kind:
                if kind == "yaw_cut" and length_columns:
                    values = [row[c] for c in length_columns if c in row and c != col]
                else:
                    values = [v for c, v in row.items() if c > col and v is not None and
                              (not isinstance(v, str) or re.fullmatch(r"\s*[-+.\d]+(?:\s*[a-zµμ°]+)?\s*", v))]
                if len(values) == 1:
                    result.append((kind, values[0], f"{obj.Label}.{col}{number} ({caption})"))
    return result


def read_native_data(master):
    """Capture scalar properties, explicit key/value rows and released cut claims."""
    data = []
    for doc in master.getDependentDocuments(True):
        for obj in doc.Objects:
            fields, cells, labels = {}, {}, {}
            if obj.isDerivedFrom("Spreadsheet::Sheet"):
                for prop in obj.getNonEmptyCells():
                    value = _scalar(getattr(obj, prop, None))
                    if value is None:
                        continue
                    cells[prop] = value
                    alias = obj.getAlias(prop)
                    if alias:
                        fields[alias] = value
                # The supplied input table has one key and one value on a
                # row. Recognize that layout even when the internal alias is
                # renamed. A heading on a different row is never a value.
                for cell, caption in cells.items():
                    m = re.fullmatch(r"A([1-9]\d*)", cell)
                    if m and isinstance(caption, str) and "B" + m[1] in cells:
                        labels[caption] = cells["B" + m[1]]
            else:
                for prop in obj.PropertiesList:
                    if prop in {"Label", "Label2", "Visibility", "ExpressionEngine"}:
                        continue
                    value = _scalar(getattr(obj, prop, None))
                    if value is not None:
                        fields[prop] = value
            if fields or cells:
                data.append({"names": [obj.Name, obj.Label], "master": doc == master,
                             "fields": fields, "labels": labels, "shop": _published_rows(obj, cells)})
    return data


def magnitude(value, expected_kind=None):
    if isinstance(value, dict):
        kind = value["unit"]
        if expected_kind is not None and kind not in ("1", expected_kind):
            raise ValueError(f"Expected {expected_kind}, found {kind}")
        return value["value"]
    return value


def unit_for(name):
    if name in {"TurnsPerCoil", "CoilType", "VaneClearanceFraction"} or name.startswith("Number"):
        return "1"
    if "Angle" in name:
        return "Angle"
    if "Weight" in name:
        return "Mass"
    if any(w in name for w in ("Diameter", "Radius", "Length", "Width", "Thickness", "Offset",
                               "Height", "Clearance", "Distance", "Margin", "MetalLength", "MetalThickness")):
        return "Length"
    return None


class NativeValues:
    def __init__(self, data, original=None, scope=None):
        self.data = data
        self.original = original
        self.scope = scope

    def field(self, name):
        matches = defaultdict(list)
        for obj in self.data:
            local = self.scope in map(_key, obj["names"])
            for channel in ("fields", "labels"):
                if channel == "labels" and not local:
                    continue
                for prop, value in obj[channel].items():
                    if _key(prop) == _key(name):
                        matches[local, channel == "fields", obj["master"], prop == name].append(value)
        if not matches:
            raise KeyError(name)
        found = matches[max(matches)]
        if not all(v == found[0] for v in found):
            raise ValueError(f"Conflicting native values for {self.scope}.{name}")
        return found[0]

    def __getattr__(self, name):
        if self.scope is None:
            return NativeValues(self.data, self.original, _key(name))
        try:
            value = self.field(name)
        except KeyError:
            # A release need not keep private field names. Evaluate the public
            # schedule and material separately and substitute only its observed
            # native fields when replaying the comparator schedule.
            if self.original is None:
                raise ValueError(f"Unresolved comparator field {self.scope}.{name}")
            value = NativeValues(self.original, scope=self.scope).field(name)
        return magnitude(value, unit_for(name))


def native_differences(data, reference, source, fraction=0.10):
    """Reconcile the stated preset, native dimensions and explicit cut directions."""
    import sys
    from unittest.mock import patch
    from comparator import _reading, table_differences
    sys.path.insert(0, str(source.parent))
    try:
        from openafpm_cad_core import dimension_tables
    finally:
        sys.path.pop(0)
    master = NativeValues(data, reference["native_data"])
    errors = []
    for group in ("magnafpm", "furling", "user"):
        for key, value in reference["config"][group].items():
            try:
                actual = magnitude(master.Spreadsheet.field(key), unit_for(key))
                same = actual == value
                if isinstance(actual, (float, int)) and isinstance(value, (float, int)):
                    same = same or abs(actual - value) <= max(math.ulp(actual), math.ulp(value))
                if not same:
                    errors.append(f"{key}: {actual!r} versus {value!r}")
            except KeyError:
                # The preset is public but its native aliases are not prescribed.
                continue
            except ValueError as error:
                errors.append(f"{key}: {error}")

    tables = reference["tables"]
    nominal = NativeValues(reference["native_data"])
    expected_cut = nominal.HighEndStop.YawPipeLength
    try:
        for obj in data:
            claims = list(obj["shop"])
            claims.extend(("fraction", v, key) for key, v in obj["fields"].items()
                          if _key(key) == "vaneclearancefraction")
            for kind, field, origin in claims:
                if isinstance(field, str):
                    readings = _reading(field, "mm" if kind == "yaw_cut" else None)
                    units = "length" if kind == "yaw_cut" else "number"
                    if set(readings) != {units} or len(readings[units]) != 1:
                        errors.append(f"Incompatible value for native {kind}: {field}")
                        continue
                    field = float(readings[units][0][0])
                if kind == "fraction":
                    if magnitude(field, "1") != fraction:
                        errors.append("Native clearance fraction disagrees with the requested build")
                elif abs(round(expected_cut) - magnitude(field, "Length")) > 0.5:
                    errors.append(f"{origin} gives a native yaw-pipe cut of "
                                  f"{magnitude(field, 'Length')} mm, expected {round(expected_cut)} mm displayed")
        # Resin consumption additionally needs the alternator geometry; its
        # shop values are checked directly through the independent public call.
        with patch.object(dimension_tables, "create_resin_table", lambda *_: tables[-1]):
            measured = dimension_tables.get_dimension_tables(master, None)
        errors.extend(table_differences(tables, measured))
    except (ValueError, TypeError, KeyError, AttributeError, ArithmeticError) as error:
        errors.append(f"Cannot reconcile the native manufacturing schedule: {error}")
    return errors
