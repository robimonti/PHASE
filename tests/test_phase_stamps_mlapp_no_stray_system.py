import re
import zipfile
from pathlib import Path

def _read_xml(path: Path) -> str:
    with zipfile.ZipFile(path, "r") as z:
        return z.read("matlab/document.xml").decode("utf-8")

def test_every_system_call_is_inside_isunix_block(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    for m in re.finditer(r"system\(", xml):
        start = max(0, m.start() - 100)
        preceding = xml[start:m.start()]
        assert "if isunix" in preceding or "else" in preceding, (
            f"system() at offset {m.start()} not inside an isunix branch: "
            f"...{preceding[-60:]}system()..."
        )
