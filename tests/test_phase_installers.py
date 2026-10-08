"""Checks for the PHASE 7 per-user installers without installing software."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def load_unix_installer(phase_root: Path):
    path = phase_root / "installer" / "install-phase-unix.py"
    spec = importlib.util.spec_from_file_location("phase_unix_installer", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_unix_installer_copies_only_runtime_files(phase_root, tmp_path):
    installer = load_unix_installer(phase_root)
    installer.SOURCE_FOR_COPY = phase_root
    ignored = installer.ignored(str(phase_root), [
        ".git", "tests", "installer", "PHASE_Hub.m", "notes.docx"
    ])
    assert ignored == {".git", "tests", "installer", "notes.docx"}
    assert installer.ignored(str(phase_root / "tools"), [
        "srtm_precache.py", "snap_dim_version_check.py", "dev.py"
    ]) == {"dev.py"}
    installer.validate_runtime(phase_root)


def test_unix_installer_preserves_unmanaged_destination(phase_root, tmp_path):
    installer = load_unix_installer(phase_root)
    prefix = tmp_path / "existing"
    prefix.mkdir()
    (prefix / "my-file.txt").write_text("keep me", encoding="utf-8")
    with pytest.raises(RuntimeError, match="not marked as a managed"):
        installer.ensure_safe_prefix(prefix, phase_root)
    assert (prefix / "my-file.txt").read_text(encoding="utf-8") == "keep me"


def test_unix_launcher_targets_installed_hub(phase_root, tmp_path):
    installer = load_unix_installer(phase_root)
    launcher = tmp_path / "launch-phase.sh"
    prefix = tmp_path / "PHASE with spaces"
    installer.write_launcher(
        launcher, prefix, Path("/opt/MATLAB/bin/matlab"),
        Path("/opt/esa-snap/bin/gpt"), Path("/opt/python/bin/python")
    )
    content = launcher.read_text(encoding="utf-8")
    assert "PHASE_Hub" in content
    assert "PHASE_GPTBIN=" in content
    assert "PHASE_PYTHON=" in content
    assert "PHASE with spaces/engine" in content
    assert launcher.stat().st_mode & 0o111


def test_windows_installer_has_one_hub_shortcut(phase_root):
    script = (phase_root / "installer" / "install-phase.ps1").read_text(encoding="utf-8-sig")
    assert "@{ Name = 'PHASE'; Launcher = 'PHASE_Hub.m'; Function = 'PHASE_Hub' }" in script
    assert "@{ Name = 'PHASE Preprocessing'" not in script
    assert "@{ Name = 'PHASE StaMPS'" not in script
    assert "@{ Name = 'PHASE Model'" not in script
    assert "Set-PhaseGptEnvVar -SnapGpt" in script
    assert "Join-Path $desktop 'PHASE.lnk'" in script
    assert "Join-Path $desktop 'PHASE 7.lnk'" in script
    assert "Removed old PHASE 7 desktop shortcut" in script
