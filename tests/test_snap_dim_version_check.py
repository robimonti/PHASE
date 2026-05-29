"""Tests for tools/snap_dim_version_check.py — offline."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import snap_dim_version_check as sdv


# ---- detect_dim_snap_version --------------------------------------------

def test_detect_dim_version_module_version_snap9(tmp_path: Path):
    """SNAP writes the producing version as `moduleVersion` on every
    Processing_Graph node — this is what real `.dim` files contain."""
    f = tmp_path / "p.dim"
    f.write_text(
        '<Dimap_Document>\n'
        '  <Dataset_Sources><MDElem name="metadata">\n'
        '    <MDElem name="Processing_Graph">\n'
        '      <MDElem name="node.0">\n'
        '        <MDATTR name="moduleVersion" type="ascii" mode="rw">9.0.0</MDATTR>\n'
        '      </MDElem>\n'
        '      <MDElem name="node.1">\n'
        '        <MDATTR name="moduleVersion" type="ascii" mode="rw">9.0.0</MDATTR>\n'
        '      </MDElem>\n'
        '    </MDElem>\n'
        '  </MDElem></Dataset_Sources>\n'
        '</Dimap_Document>\n',
        encoding="utf-8",
    )
    assert sdv.detect_dim_snap_version(f) == "9.0.0"


def test_detect_dim_version_module_version_snap13(tmp_path: Path):
    f = tmp_path / "p.dim"
    f.write_text(
        '<Dimap_Document>\n'
        '  <MDElem name="Processing_Graph">\n'
        '    <MDElem name="node.0">\n'
        '      <MDATTR name="moduleVersion" type="ascii" mode="rw">13.0.0</MDATTR>\n'
        '    </MDElem>\n'
        '  </MDElem>\n'
        '</Dimap_Document>\n',
        encoding="utf-8",
    )
    assert sdv.detect_dim_snap_version(f) == "13.0.0"


def test_detect_dim_version_module_version_picks_highest(tmp_path: Path):
    """Robustness: if a third-party operator embeds an older moduleVersion,
    the SNAP-framework nodes (highest version) win."""
    f = tmp_path / "p.dim"
    f.write_text(
        '<MDElem name="Processing_Graph">\n'
        '  <MDATTR name="moduleVersion" type="ascii">13.0.0</MDATTR>\n'
        '  <MDATTR name="moduleVersion" type="ascii">2.5.1</MDATTR>\n'
        '  <MDATTR name="moduleVersion" type="ascii">13.0.0</MDATTR>\n'
        '</MDElem>\n',
        encoding="utf-8",
    )
    assert sdv.detect_dim_snap_version(f) == "13.0.0"


def test_detect_dim_version_real_snap9_artifact():
    """Smoke-test against an actual SNAP 9 .dim from the Bolsena dataset."""
    real = (Path(__file__).resolve().parent.parent
            / "PHASE_fork" / "PHASE_Preprocessing" / "master"
            / "S1A_IW_SLC__1SDV_20240715T052011_20240715T052038_054767_06AB2A_6A8C_split_IW1_Orb.dim")
    if not real.is_file():
        pytest.skip(f"real artifact not present: {real}")
    assert sdv.detect_dim_snap_version(real) == "9.0.0"


def test_detect_dim_version_real_snap13_artifact():
    """Smoke-test against an actual SNAP 13 .dim from a previous run."""
    real = (Path(__file__).resolve().parent.parent
            / "test_snap13" / "Phase" / "PHASE_fork" / "PHASE_Preprocessing"
            / "coreg" / "_20240703.dim")
    if not real.is_file():
        pytest.skip(f"real artifact not present: {real}")
    assert sdv.detect_dim_snap_version(real) == "13.0.0"


def test_detect_dim_version_processing_software(tmp_path: Path):
    f = tmp_path / "p.dim"
    f.write_text(
        '<Dimap_Document>\n'
        '  <Processing_Graph>\n'
        '    <node>\n'
        '      <MDATTR name="processing_software_name" type="ascii">SNAP</MDATTR>\n'
        '      <MDATTR name="processing_software_version" type="ascii">9.0.0</MDATTR>\n'
        '    </node>\n'
        '  </Processing_Graph>\n'
        '</Dimap_Document>\n',
        encoding="utf-8",
    )
    assert sdv.detect_dim_snap_version(f) == "9.0.0"


def test_detect_dim_version_snap_version_alt(tmp_path: Path):
    f = tmp_path / "p.dim"
    f.write_text(
        '<Dimap_Document>\n'
        '  <MDATTR name="snap_version" type="ascii">13.0.0</MDATTR>\n'
        '</Dimap_Document>\n',
        encoding="utf-8",
    )
    assert sdv.detect_dim_snap_version(f) == "13.0.0"


def test_detect_dim_version_missing(tmp_path: Path):
    f = tmp_path / "p.dim"
    f.write_text('<Dimap_Document><foo/></Dimap_Document>\n', encoding="utf-8")
    assert sdv.detect_dim_snap_version(f) is None


def test_detect_dim_version_nonexistent(tmp_path: Path):
    assert sdv.detect_dim_snap_version(tmp_path / "missing.dim") is None


# ---- major --------------------------------------------------------------

def test_major_extraction():
    assert sdv.major("9.0.0") == 9
    assert sdv.major("13.0.0") == 13
    assert sdv.major("13.0") == 13
    assert sdv.major("13") == 13


# ---- check_compatibility (with injected paths/mocks) --------------------

def _write_dim(p: Path, version: str) -> None:
    p.write_text(
        f'<MDATTR name="processing_software_version" type="ascii">{version}'
        '</MDATTR>\n',
        encoding="utf-8",
    )


def test_check_compatibility_all_match(tmp_path: Path, monkeypatch):
    d1 = tmp_path / "a.dim"; _write_dim(d1, "13.0.0")
    d2 = tmp_path / "b.dim"; _write_dim(d2, "13.0.0")
    monkeypatch.setattr(sdv, "detect_gpt_snap_version", lambda _p: "13.0.0")
    report = sdv.check_compatibility([d1, d2], tmp_path / "gpt.exe")
    assert len(report["matched"]) == 2
    assert report["mismatched"] == []
    assert report["unknown"] == []
    assert report["gpt_version"] == "13.0.0"


def test_check_compatibility_mismatch(tmp_path: Path, monkeypatch):
    d_old = tmp_path / "old.dim"; _write_dim(d_old, "9.0.0")
    d_new = tmp_path / "new.dim"; _write_dim(d_new, "13.0.0")
    monkeypatch.setattr(sdv, "detect_gpt_snap_version", lambda _p: "13.0.0")
    report = sdv.check_compatibility([d_old, d_new], tmp_path / "gpt.exe")
    assert len(report["matched"]) == 1
    assert len(report["mismatched"]) == 1
    mism_path, dim_ver, gpt_ver = report["mismatched"][0]
    assert "old.dim" in mism_path
    assert dim_ver == "9.0.0"
    assert gpt_ver == "13.0.0"


def test_check_compatibility_unknown_version(tmp_path: Path, monkeypatch):
    d = tmp_path / "novers.dim"
    d.write_text("<no version info>", encoding="utf-8")
    monkeypatch.setattr(sdv, "detect_gpt_snap_version", lambda _p: "13.0.0")
    report = sdv.check_compatibility([d], tmp_path / "gpt.exe")
    assert len(report["unknown"]) == 1


def test_check_compatibility_skips_when_gpt_unknown(tmp_path: Path,
                                                    monkeypatch):
    d = tmp_path / "a.dim"; _write_dim(d, "9.0.0")
    monkeypatch.setattr(sdv, "detect_gpt_snap_version", lambda _p: None)
    report = sdv.check_compatibility([d], tmp_path / "gpt.exe")
    # No check performed: skip silently rather than break user workflows.
    assert report["matched"] == []
    assert report["mismatched"] == []
    assert report["unknown"] == []
    assert report["gpt_version"] is None


# ---- warn_if_mismatch ---------------------------------------------------

def test_warn_if_mismatch_silent_when_no_mismatch():
    captured: list[str] = []
    report = {"mismatched": [], "matched": [], "unknown": [],
              "gpt_version": "13.0.0"}
    fired = sdv.warn_if_mismatch(report, write=captured.append)
    assert fired is False
    assert captured == []


def test_warn_if_mismatch_prints_actionable_diagnostic():
    captured: list[str] = []
    report = {
        "mismatched": [("F:/path/old.dim", "9.0.0", "13.0.0")],
        "matched": [], "unknown": [], "gpt_version": "13.0.0",
    }
    fired = sdv.warn_if_mismatch(report, write=captured.append)
    assert fired is True
    blob = "\n".join(captured)
    assert "9.0.0" in blob
    assert "13.0.0" in blob
    assert "old.dim" in blob
    assert "re-run" in blob.lower() or "re-process" in blob.lower()
