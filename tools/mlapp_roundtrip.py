"""Unpack, edit, and repack a MATLAB .mlapp while preserving Office Open XML file order."""
import re
import shutil
import sys
import zipfile
from pathlib import Path


# Canonical file order for .mlapp (matches MATLAB expectations).
_ORDER_PREFIX = [
    "[Content_Types].xml",
    "_rels/.rels",
]


def _file_sort_key(name: str) -> tuple[int, int, str]:
    """Earlier tuple => earlier in zip."""
    if name == "[Content_Types].xml":
        return (0, 0, name)
    if name == "_rels/.rels":
        return (0, 1, name)
    if name.startswith("metadata/"):
        return (1, 0, name)
    if name == "matlab/document.xml":
        return (2, 0, name)
    if name.startswith("matlab/"):
        return (2, 1, name)
    if name.startswith("appdesigner/"):
        return (3, 0, name)
    return (9, 0, name)


def edit_mlapp(mlapp_path: Path, edit_fn) -> None:
    """Edit the document.xml of an .mlapp via edit_fn(xml_str) -> xml_str."""
    with zipfile.ZipFile(mlapp_path, "r") as zin:
        names = zin.namelist()
        contents = {n: zin.read(n) for n in names}

    xml = contents["matlab/document.xml"].decode("utf-8")
    new_xml = edit_fn(xml)
    contents["matlab/document.xml"] = new_xml.encode("utf-8")

    # Write a new .mlapp with canonical ordering.
    tmp = mlapp_path.with_suffix(".mlapp.tmp")
    ordered = sorted(contents.keys(), key=_file_sort_key)
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as zout:
        for name in ordered:
            zout.writestr(name, contents[name])
    shutil.move(tmp, mlapp_path)


if __name__ == "__main__":
    # CLI: python mlapp_roundtrip.py <file> <sed-like-pattern> <replacement>
    # Used by workers for simple regex-substitution edits.
    path, pattern, repl = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    edit_mlapp(path, lambda x: re.sub(pattern, repl, x, count=1))
