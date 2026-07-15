import re
import zipfile
from pathlib import Path

def _read_xml(path: Path) -> str:
    with zipfile.ZipFile(path, "r") as z:
        return z.read("matlab/document.xml").decode("utf-8")

def test_every_system_call_is_inside_isunix_block(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    code = re.sub(r"(?m)^\s*%.*(?:\n|$)", "\n", xml)
    # Match the MATLAB system() function, not names such as sp_system()
    # mentioned in comments or helper calls.
    for m in re.finditer(r"(?<![A-Za-z0-9_])system\(", code):
        start = max(0, m.start() - 100)
        preceding = code[start:m.start()]
        assert "if isunix" in preceding or "else" in preceding, (
            f"system() at offset {m.start()} not inside an isunix branch: "
            f"...{preceding[-60:]}system()..."
        )
