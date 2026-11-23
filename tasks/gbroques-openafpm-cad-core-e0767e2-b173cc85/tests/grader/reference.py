"""Measure the supplied release with only the two specified stop settings replaced."""
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path.cwd()))
import FreeCAD as App
from comparator import build
from solids import snapshot
from motion import _frames, _measured_members, _native_roles, _turn
from native import motion_trial
from native_fields import read_native_data


def memberships(wind, inventory, folder, cache):
    """Cache measured fixed/moving membership by installed occurrence and solid ordinal."""
    roles = _native_roles(inventory.frames)
    paths = {(*row["path"], ordinal): index for row in inventory.frames if row["leaf"]
             for ordinal, index in enumerate(row["members"])}
    states = {}
    if cache.is_file():
        states = {tuple(path): pair for path, pair in json.loads(cache.read_text())}
    if states.keys() != paths.keys():
        motion_trial(wind, inventory, roles["tail"], roles["hinge"], folder)
        (folder / "frames.json").write_text(json.dumps(inventory.frames))
        rows = _frames(folder)
        turn = _turn(rows[roles["hinge"]]["effective"])
        members = _measured_members(rows, turn, inventory.solids, range(len(inventory.solids)), folder)
        if any(fixed == moving for fixed, moving in members):
            raise ValueError("Unresolved moving/fixed membership in the comparator")
        states = {path: members[index] for path, index in paths.items()}
        cache.write_text(json.dumps(list(states.items())))
    ordered = [None] * len(inventory.solids)
    for path, index in paths.items():
        ordered[index] = states[path]
    inventory.frames[roles["tail"]]["mobility"] = ordered


def main():
    preset, angle, width, cut, filename = sys.argv[1:]
    wind, master, tables, config, nominal_boom = build(preset, float(angle), float(width), int(cut) if cut else None)
    folder = Path("record")
    folder.mkdir()
    native_data = read_native_data(master)
    inventory = snapshot(wind, folder)
    memberships(wind, inventory, folder, Path(filename))
    (folder / "frames.json").write_text(json.dumps(inventory.frames))
    with (folder / "result.json").open("w") as f:
        json.dump(dict(tables=tables, config=config, native_data=native_data, nominal_boom=nominal_boom,
                       boom_change=config["furling"]["BoomLength"] - nominal_boom), f)
    for doc in list(App.listDocuments(True).values())[::-1]:
        App.closeDocument(doc.Name)


if __name__ == "__main__":
    main()
