import re
import zipfile


def _stable_xml(phase_root):
    path = phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"
    with zipfile.ZipFile(path) as app:
        return app.read("matlab/document.xml").decode("utf-8")


def _text(path):
    return path.read_text(encoding="utf-8")


def test_stable_weed_standard_dev_uses_range_and_one_decimal(phase_root):
    xml = _stable_xml(phase_root)
    declarations = re.findall(
        r"app\.weed_standard_devEditField\.Limits\s*=\s*\[([^]]+)\];",
        xml,
    )
    assert declarations == ["0 100"]
    assert "app.weed_standard_devEditField.RoundFractionalValues = 'off';" in xml
    assert "app.weed_standard_devEditField.ValueDisplayFormat = '%.1f';" in xml


def test_beta_weed_standard_dev_matches_stable_range_and_precision(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    schema = _text(package / "schema.m")
    conversion = _text(package / "uiToConfig.m")
    validation = _text(package / "validateConfig.m")
    self_test = _text(package / "selfTest.m")
    ui = _text(
        phase_root / "PHASE_Preprocessing" / "phase_stamps_beta_ui" / "app.js"
    )

    assert "'Range 0–100; one decimal place.'" in schema
    assert "'weed_standard_dev'" in conversion
    assert "cfg.weed_standard_dev < 0 || cfg.weed_standard_dev > 100" in validation
    assert "cfg.weed_standard_dev * 10" in validation
    assert "cfg.weed_standard_dev = 37.5" in self_test
    assert 'item.id === "weed_standard_dev"' in ui
    assert 'control.max = "100"' in ui
    assert 'control.step = "0.1"' in ui
