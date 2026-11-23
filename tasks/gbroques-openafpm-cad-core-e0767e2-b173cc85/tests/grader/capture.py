"""Invoke the public entry point, measure the returned model and transfer its native graph."""
import json
import os
from pathlib import Path
import shutil
import stat
import sys


def mirror_directory(source, destination, user_files=False):
    """Copy the work tree as its owner, keeping regular project files and their layout."""
    def copy(directory, target, ancestry):
        target.mkdir(parents=True, exist_ok=True)
        try:
            entries = list(directory.iterdir())
        except PermissionError:
            return
        for entry in entries:
            try:
                path = entry.resolve(strict=True)
                info = path.stat()
                if user_files and directory == source and info.st_uid != os.getuid():
                    continue
                mode = info.st_mode
            except (OSError, RuntimeError):
                continue
            dest = target / entry.name
            if stat.S_ISREG(mode):
                try:
                    shutil.copyfile(path, dest)
                except PermissionError:
                    continue
                dest.chmod(stat.S_IMODE(mode) | stat.S_IRUSR | stat.S_IWUSR)
            elif stat.S_ISDIR(mode) and path not in ancestry and path != destination:
                copy(path, dest, ancestry | {path})
    copy(source, destination, {source})


def transfer(wind, master, work, mirror):
    mirror_directory(work, mirror)
    locations = {}
    originals = {}
    for doc in wind.getDependentDocuments(True):
        old = Path(doc.FileName).resolve() if doc.FileName else None
        originals[doc] = old
        if old is None:
            new = mirror / ".documents" / (doc.Name + ".FCStd")
        elif old.is_relative_to(work):
            new = mirror / old.relative_to(work)
        else:
            # Preserve supporting files beside a component generated in an
            # invocation-local directory outside its default temporary root.
            directory = mirror / ".documents" / old.parent.relative_to("/")
            if not directory.exists():
                mirror_directory(old.parent, directory, user_files=old.parent.stat().st_uid != os.getuid())
            new = directory / old.name
        new.parent.mkdir(parents=True, exist_ok=True)
        locations[doc] = new
    # saveCopy retains filenames and the links between them in memory. Unlike
    # saveAs on an already named document, it also preserves the saved relative
    # link layout when the complete directory is transferred together.
    for doc, new in locations.items():
        if originals[doc] is None:
            doc.saveAs(str(new))
    for doc, new in locations.items():
        if originals[doc] is not None:
            doc.saveCopy(str(new))
    if any(old is None or not old.is_relative_to(work) for old in originals.values()):
        from xml.etree import ElementTree
        from zipfile import ZipFile
        lookup = {(old or locations[doc]): locations[doc] for doc, old in originals.items()}
        lookup.update({path: path for path in locations.values()})
        for doc, path in locations.items():
            with ZipFile(path) as saved:
                tree = ElementTree.fromstring(saved.read("Document.xml"))
                changed = False
                for link in tree.iter("XLink"):
                    name = link.get("file")
                    if not name:
                        continue
                    target = lookup.get(((originals[doc] or path).parent / name).resolve())
                    if target is not None and (path.parent / name).resolve() != target:
                        link.set("file", os.path.relpath(target, path.parent))
                        changed = True
                if not changed:
                    continue
                temp = path.with_suffix(".tmp")
                with ZipFile(temp, "w") as rebuilt:
                    for item in saved.infolist():
                        content = (ElementTree.tostring(tree, encoding="utf-8") if
                                   item.filename == "Document.xml" else saved.read(item))
                        rebuilt.writestr(item, content)
            temp.replace(path)
    return {key: str(locations[doc].relative_to(mirror)) for key, doc in
            (("document", wind), ("master", master))}


def invoke(mode, preset, fraction, cut, reader, output):
    from openafpm_cad_core.get_default_parameters import get_default_parameters
    from openafpm_cad_core.load import load_turbine
    from openafpm_cad_core.furl_transform import get_furl_transform
    from openafpm_cad_core.dimension_tables import load_dimension_tables

    output = Path(output)
    cfg = get_default_parameters(preset)
    args = [cfg["magnafpm"], dict(cfg["furling"]), cfg["user"]]
    result = {"preset": json.loads(json.dumps(args))}
    if cut:
        args[1]["BoomLength"] = int(cut)
    if fraction != "default":
        args[1]["VaneClearanceFraction"] = float(fraction)
    if mode == "tables":
        result["tables"] = load_dimension_tables(*args)
    else:
        wind, master = load_turbine(*args)
        import FreeCAD as App
        if not all(isinstance(doc, App.Document) for doc in (wind, master)):
            raise ValueError("The assembly call must return two native documents")
        if master == wind or master not in wind.getDependentDocuments(True):
            raise ValueError("The returned master is outside the production document dependencies")
        sys.path.insert(0, reader)
        from native import governing_settings
        from solids import snapshot
        snapshot(wind, output / "returned")
        _, result["returned_angle"], result["returned_width"] = governing_settings(wind, master)
        result["motion"] = get_furl_transform(wind, master)
        result.update(transfer(wind, master, Path.cwd(), output / "mirror"))
    (output / "result.json").write_text(json.dumps(result, allow_nan=False))
