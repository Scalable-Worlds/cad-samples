"""Snapshot the complete installed material and its independently mounted frames."""
import json
import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import FreeCAD
import Part
from frames import Inventory


def snapshot(doc, directory):
    folder = Path(directory)
    folder.mkdir(exist_ok=True, parents=True)
    inventory = Inventory(doc)
    if not inventory.solids:
        raise ValueError("The production document has no installed closed solids")
    Part.makeCompound(inventory.solids).exportBrep(str(folder / "assembly.brep"))
    (folder / "datums.json").write_text(json.dumps(dict(quaternions=inventory.datums, invalid=inventory.invalid)))
    (folder / "frames.json").write_text(json.dumps(inventory.frames))
    return inventory
