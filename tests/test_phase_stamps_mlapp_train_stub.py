import zipfile
from pathlib import Path

def _read_xml(path):
    with zipfile.ZipFile(path) as z:
        return z.read("matlab/document.xml").decode("utf-8")

def test_train_has_windows_empty_branch(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "APS_CONFIG.sh" in xml
    assert "train_path = ''" in xml


def test_canonical_warning_id_present(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "StaMPS:phase:trainNotAvailable" in xml
