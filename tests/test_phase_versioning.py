"""The hub and module sidebars must display the same release version."""

import re


def test_all_module_versions_use_shared_release_metadata(phase_root):
    resolver = (phase_root / "+phase_project" / "displayVersion.m").read_text()
    assert "install.json" in resolver
    assert "phase-release.json" in resolver
    sources = (
        phase_root / "+phase_hub" / "App.m",
        phase_root / "PHASE_Preprocessing" / "+phase_preprocessing_beta" / "App.m",
        phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta" / "App.m",
        phase_root / "+phase_model_beta" / "App.m",
    )
    for source in sources:
        text = source.read_text(encoding="utf-8")
        assert "phase_project.displayVersion(" in text
        assert "'6.1.0'" not in text
    updater = (phase_root / "phase_update.py").read_text(encoding="utf-8")
    assert '"+phase_project/displayVersion.m"' in updater


def test_development_version_matches_installer_previews(phase_root):
    resolver = (phase_root / "+phase_project" / "displayVersion.m").read_text()
    match = re.search(r"value = 'v(\d+\.\d+\.\d+) dev';", resolver)
    assert match
    version = match.group(1)
    assert f"v{version} preview" in (
        phase_root / "installer" / "PHASEInstaller.swift"
    ).read_text()
    assert f"v{version} preview" in (
        phase_root / "installer" / "install-phase.ps1"
    ).read_text()
    assert f"version    = '{version}.0'" in (
        phase_root / "installer" / "compile-to-exe.ps1"
    ).read_text()
