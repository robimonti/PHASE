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
    assert "No legacy MLAPP is copied" in source
    assert "PHASE StaMPS is opened by preprocessing" in source


def test_beta_installer_clones_branch_and_launches_standalone_m_files(phase_root):
    source = _installer(phase_root)

    assert "[string]$PhaseBranch = 'codex/phase-stamps-beta'" in source
    assert "$Script:PhaseBranch = $PhaseBranch" in source
    assert "PHASE_Preprocessing_beta.m" in source
    assert "PHASE_Model_beta.m" in source
    assert "$sc.TargetPath = $MatlabExe" in source
    assert "$($a.Function)" in source
    assert "Open $($a.Name) in MATLAB App Designer" not in source


def test_beta_installer_cleans_legacy_runtime_and_hides_engine(phase_root):
    source = _installer(phase_root)

    assert "function Remove-PhaseLegacyRuntimeFiles" in source
    assert "'PHASE_Preprocessing.mlapp'" in source
    assert "'PHASE_Preprocessing\\PHASE_StaMPS.mlapp'" in source
    assert "'PHASE_model.mlapp'" in source
    assert "'srtm_precache.py', 'snap_dim_version_check.py'" in source
    assert "$_.Name -notin $runtimeTools" in source
    assert "Retained SNAP runtime helpers in tools." in source
    assert "function Set-PhaseEngineHidden" in source
    assert "[System.IO.FileAttributes]::Hidden" in source
    assert "Set-PhaseEngineHidden -EngineDir $engineDir" in source
    assert "function Assert-PhaseStandaloneRuntime" in source
    assert "Assert-PhaseStandaloneRuntime -PhaseDir $phaseDir" in source
    assert "Nothing was cleaned" in source


def test_installer_compiler_resolves_defaults_after_parameter_binding(phase_root):
    compiler = (phase_root / "installer" / "compile-to-exe.ps1").read_text(encoding="utf-8")

    param_block = compiler[compiler.index("param(") : compiler.index(")\n\n$ErrorActionPreference")]
    assert "Join-Path $PSScriptRoot" not in param_block
    assert "$scriptPath = $MyInvocation.MyCommand.Path" in compiler
    assert "$Source = Join-Path $scriptDir 'install-phase.ps1'" in compiler
    assert "$Output = Join-Path $scriptDir 'install-phase-beta.exe'" in compiler
