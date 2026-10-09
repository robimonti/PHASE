"""Checks for the PHASE 7 per-user installers without installing software."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil
import subprocess

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
    assert "export STAMPS=" in content
    assert "export APS_toolbox=" in content
    assert "/opt/MATLAB/bin:" in content
    assert "StaMPS/matlab" in content
    assert "TRAIN/matlab" in content
    assert launcher.stat().st_mode & 0o111


def test_unix_runtime_configs_do_not_use_upstream_example_paths(phase_root, tmp_path):
    installer = load_unix_installer(phase_root)
    stage = tmp_path / "stage"
    (stage / "engine" / "StaMPS").mkdir(parents=True)
    (stage / "engine" / "TRAIN").mkdir(parents=True)
    prefix = tmp_path / "PHASE with spaces"
    installer.configure_unix_runtimes(stage, prefix)
    stamps = (stage / "engine" / "StaMPS" / "StaMPS_CONFIG.bash").read_text()
    train = (stage / "engine" / "TRAIN" / "APS_CONFIG.sh").read_text()
    assert str(prefix / "engine" / "StaMPS") in stamps
    assert str(prefix / "engine" / "TRAIN") in train
    assert "/home/ahooper" not in stamps
    assert "/nfs/see-fs" not in train
    if shutil.which("bash"):
        subprocess.run(["bash", "-n", str(stage / "engine" / "StaMPS" / "StaMPS_CONFIG.bash")], check=True)
        subprocess.run(["bash", "-n", str(stage / "engine" / "TRAIN" / "APS_CONFIG.sh")], check=True)


def test_macos_runtime_builder_never_writes_inside_source(phase_root, tmp_path, monkeypatch):
    path = phase_root / "installer" / "prepare-macos-runtime.py"
    spec = importlib.util.spec_from_file_location("phase_macos_runtime", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    monkeypatch.setattr(module.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(module.platform, "machine", lambda: "arm64")
    stamps = tmp_path / "StaMPS"
    train = tmp_path / "TRAIN"
    stamps.mkdir()
    train.mkdir()
    with pytest.raises(RuntimeError, match="must be separate"):
        module.prepare(stamps, train, tmp_path / "snaphu", tmp_path / "triangle",
                       tmp_path / "gawk",
                       stamps / "runtime")
    assert not (stamps / "runtime").exists()


def test_macos_installer_reports_missing_native_psi_tools(phase_root, tmp_path):
    installer = load_unix_installer(phase_root)
    stamps = tmp_path / "StaMPS"
    stamps.mkdir()
    missing = installer.macos_psi_missing(stamps, None)
    assert "StaMPS/bin/calamp" in missing
    assert "snaphu" in missing
    assert "triangle" in missing
    assert "gawk" in missing
    assert any("TRAIN runtime" in item for item in missing)


def test_macos_installer_accepts_complete_arm64_psi_runtime(phase_root, tmp_path, monkeypatch):
    installer = load_unix_installer(phase_root)
    stamps = tmp_path / "StaMPS"
    train = tmp_path / "TRAIN"
    for name in ("calamp", "cpxsum", "pscphase", "pscdem", "psclonlat",
                 "selpsc_patch", "selsbc_patch"):
        binary = stamps / "bin" / name
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_bytes(b"arm64 binary")
        binary.chmod(0o755)
    for name in ("snaphu", "triangle", "gawk"):
        binary = stamps / "external" / name / "bin" / name
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_bytes(b"arm64 binary")
        binary.chmod(0o755)
    (train / "matlab").mkdir(parents=True)
    (train / "matlab" / "aps_linear.m").write_text("", encoding="utf-8")
    monkeypatch.setattr(installer.subprocess, "run", lambda *args, **kwargs:
                        subprocess.CompletedProcess(args[0], 0, "arm64\n", ""))
    assert installer.macos_psi_missing(stamps, train) == []


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
    assert 'x:Name="LaunchPhaseBtn" Content="Launch PHASE"' in script
    assert "Start-Process -FilePath $launcher -ErrorAction Stop" in script
    assert "Join-Path (Join-Path $Script:State.InstallDir 'PHASE') 'PHASE.lnk'" in script
    assert "OpenFolderBtn" not in script
    assert "v7.0.0 preview" in script
    assert "Roberto Monti · pyccino" in script
    assert "Install once. Work on any number of projects." in script
    assert "Each project can live in any folder or drive" in script
    assert "Persistent scatterer Highly Automated Suite for Environmental monitoring" in script
    assert "$Script:EmbeddedLogoBase64 = ''" in script
    compiler = (phase_root / "installer" / "compile-to-exe.ps1").read_text(encoding="utf-8-sig")
    assert "version    = '7.0.0.0'" in compiler
    assert "company    = 'Roberto Monti and pyccino'" in compiler
    assert "[Convert]::ToBase64String([IO.File]::ReadAllBytes($logoPath))" in compiler
