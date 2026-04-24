import zipfile
from pathlib import Path


def _read_xml(path):
    with zipfile.ZipFile(path) as z:
        return z.read("matlab/document.xml").decode("utf-8")


def test_config_lookup_has_isunix_branch(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "if isunix" in xml
    assert "StaMPS_CONFIG.ps1" in xml and "StaMPS_CONFIG.bash" in xml
