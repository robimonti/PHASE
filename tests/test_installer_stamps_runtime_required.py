from pathlib import Path


def _installer(phase_root: Path) -> str:
    return (phase_root / "installer" / "install-phase.ps1").read_text(encoding="utf-8")


def test_missing_stamps_native_binaries_are_fatal(phase_root):
    source = _installer(phase_root)

    failed_branch = source.index("if ($binOk)")
    end = source.index("# Task 6: GMT", failed_branch)
    block = source[failed_branch:end]

    assert "9 mandatory StaMPS binaries ready" in block
    assert "step 6 cannot run without snaphu.exe" in block
    assert 'throw "The mandatory StaMPS Windows binaries could not be installed.' in block


def test_installer_does_not_advertise_invalid_global_stamps_shortcut(phase_root):
    source = _installer(phase_root)

    apps_start = source.index("$apps = @(")
    apps_end = source.index("foreach ($a in $apps)", apps_start)
    apps = source[apps_start:apps_end]

    assert "@{ Name = 'PHASE StaMPS'" not in apps
    assert "Removed obsolete global PHASE StaMPS shortcut" in apps
    assert "To resume a dataset later, open PHASE_StaMPS.mlapp inside" in source
    assert "PHASE StaMPS is opened by preprocessing" in source
