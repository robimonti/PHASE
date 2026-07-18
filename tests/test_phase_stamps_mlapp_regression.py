import zipfile
from pathlib import Path


def _read_xml(path: Path) -> str:
    with zipfile.ZipFile(path) as mlapp_zip:
        return mlapp_zip.read("matlab/document.xml").decode("utf-8")


def test_stale_app_designer_save_did_not_return(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "AutoDetectParameters(app);is" not in xml
    assert "% begin TRAIN interactive recovery" in xml
    assert "self-bootstrapping .bat shim" in xml


def test_train_and_subtr_tropo_must_both_request_correction(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "tropo_correction_enabled = (train_flag == 0) && ..." in xml
    assert (
        "train_requested = (train_flag == 0) && "
        "strcmpi(strtrim(subtr_tropo), 'y');" in xml
    )
