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


def test_installer_creates_dataset_safe_stamps_shortcut(phase_root):
    source = _installer(phase_root)

    apps_start = source.index("$apps = @(")
    apps_end = source.index("foreach ($a in $apps)", apps_start)
    apps = source[apps_start:apps_end]

    assert "@{ Name = 'PHASE StaMPS'" in apps
    assert "DatasetScoped = $true" in apps
    assert "datasetDir = uigetdir" in source
    assert "$($a.Function)(datasetDir)" in source
    assert "if ~isequal(datasetDir,0)" in source
    assert "No legacy MLAPP is copied" in source
    assert '"PHASE StaMPS.lnk"' in source


def test_installer_clones_main_and_launches_production_m_files(phase_root):
    source = _installer(phase_root)

    assert "[string]$PhaseBranch = 'main'" in source
    assert "$Script:PhaseBranch = $PhaseBranch" in source
    assert "Launcher = 'PHASE_Preprocessing.m'; Function = 'PHASE_Preprocessing'" in source
    assert "Launcher = 'PHASE_Model.m'; Function = 'PHASE_Model'" in source
    assert "$sc.TargetPath = $MatlabExe" in source
    assert "$($a.Function)" in source
    assert "Open $($a.Name) in MATLAB App Designer" not in source


def test_installer_cleans_legacy_runtime_and_keeps_engine_editable(phase_root):
    source = _installer(phase_root)

    assert "function Remove-PhaseLegacyRuntimeFiles" in source
    assert "'legacy'" in source
    assert "'srtm_precache.py', 'snap_dim_version_check.py'" in source
    assert "$_.Name -notin $runtimeTools" in source
    assert "Retained SNAP runtime helpers in tools." in source
    assert "function Set-PhaseEngineHidden" not in source
    assert "Set-PhaseEngineHidden -EngineDir $engineDir" not in source
    assert "function Set-PhaseEngineVisible" in source
    assert "-band (-bnot [System.IO.FileAttributes]::Hidden)" in source
    assert "Set-PhaseEngineVisible -EngineDir $engineDir" in source
    assert "editable MATLAB sources live in the visible folder" in source
    assert "function Assert-PhaseStandaloneRuntime" in source
    assert "Assert-PhaseStandaloneRuntime -PhaseDir $phaseDir" in source
    assert "Nothing was cleaned" in source
    assert "geoSplinter\\windows\\geoSplinter_analysis.exe" in source
    assert "geoSplinter\\windows\\geoSplinter_synthesis.exe" in source


def test_installer_compiler_resolves_defaults_after_parameter_binding(phase_root):
    compiler = (phase_root / "installer" / "compile-to-exe.ps1").read_text(encoding="utf-8")

    param_block = compiler[compiler.index("param(") : compiler.index(")\n\n$ErrorActionPreference")]
    assert "Join-Path $PSScriptRoot" not in param_block
    assert "$scriptPath = $MyInvocation.MyCommand.Path" in compiler
    assert "$Source = Join-Path $scriptDir 'install-phase.ps1'" in compiler
    assert "$Output = Join-Path $scriptDir 'install-phase.exe'" in compiler


def test_production_launchers_wrap_the_validated_standalone_engines(phase_root):
    prep = (phase_root / "PHASE_Preprocessing.m").read_text(encoding="utf-8")
    stamps = (phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.m").read_text(encoding="utf-8")
    model = (phase_root / "PHASE_Model.m").read_text(encoding="utf-8")

    assert "PHASE_Preprocessing_beta()" in prep
    assert "PHASE_StaMPS_beta(workDir)" in stamps
    assert "PHASE_Model_beta()" in model
    assert not (phase_root / "PHASE_Preprocessing.mlapp").exists()
    assert not (phase_root / "PHASE_model.mlapp").exists()
    assert not (phase_root / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp").exists()
    assert (phase_root / "legacy" / "PHASE_Preprocessing.mlapp").is_file()
    assert (phase_root / "legacy" / "PHASE_model.mlapp").is_file()
    assert (phase_root / "legacy" / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp").is_file()
    notes = (phase_root / "RELEASE_NOTES.md").read_text(encoding="utf-8")
    readme = (phase_root / "README.md").read_text(encoding="utf-8")
    assert "Do not clone the repository" in notes
    assert "releases/latest/download/install-phase.exe" in readme
    assert not (phase_root / "installer" / "install-phase.exe").exists()
    assert (phase_root / "legacy" / "install-phase-pre-v6.exe").is_file()
