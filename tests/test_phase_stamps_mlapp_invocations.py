import zipfile
from pathlib import Path


def _read_xml(path):
    with zipfile.ZipFile(path) as z:
        return z.read("matlab/document.xml").decode("utf-8")


def test_all_four_system_calls_have_isunix_branch(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    # Expect >=6 `if isunix` blocks (1 from 5.1-A config + 1 from 5.1-B train + 4 from this)
    assert xml.count("if isunix") >= 6
    assert "mt_prep_snap.bat" in xml
    assert "self-bootstrapping .bat shim" in xml
