"""Identify additional native builds by their governing settings and linked documents."""
from xml.etree import ElementTree
from zipfile import BadZipFile, ZipFile, is_zipfile


def conflicting_assemblies(out, production):
    archive = {}
    for path in out.rglob("*"):
        if not path.is_file() or not is_zipfile(path):
            continue
        try:
            with ZipFile(path) as file:
                document = ElementTree.fromstring(file.read("Document.xml"))
        except (BadZipFile, KeyError, OSError, ElementTree.ParseError):
            continue
        targets = {(path.parent / link.get("file")).resolve() for link in document.iter("XLink")
                   if link.get("file")}
        governs = any((item.get("alias") or item.get("name")) in
                      ("MaximumFurlAngle", "HighEndStopWidth")
                      for item in document.iter() if item.tag in ("Cell", "Property"))
        archive[path.resolve()] = targets, governs
    def closure(path):
        pending = [path]
        seen = set()
        while pending:
            here = pending.pop()
            if here in archive and here not in seen:
                seen.add(here)
                pending.extend(archive[here][0] - seen)
        return seen
    released = closure(production.resolve())
    production_links = archive.get(production.resolve(), (set(), False))[0]
    extra = []
    for path in archive.keys() - released:
        links, governs = archive[path]
        duplicates_root = bool(production_links and production_links <= links)
        if governs or duplicates_root:
            extra.append(path.relative_to(out))
    return sorted(extra)
