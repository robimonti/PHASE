"""Regression guards for the native Windows Triangle dependencies."""


def test_installer_bundles_mingw_runtime_and_probes_triangle(phase_root):
    source = (phase_root / "installer" / "install-phase.ps1").read_text(encoding="utf-8-sig")
    for name in ("libstdc++-6.dll", "libgcc_s_seh-1.dll", "libwinpthread-1.dll"):
        assert name in source
    assert "PHASE\\portable-git\\mingw64\\bin" in source
    assert "Install-TriangleRuntime -TriangleBinDir $triangleBinDir" in source
    assert "Start-Process -FilePath $triangle -ArgumentList '-e', 'phase-probe.1.node'" in source
    assert "phase-probe.2.edge" in source


def test_managed_install_repairs_triangle_before_stamps_runs(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    runtime = (package / "prepareRuntime.m").read_text(encoding="utf-8")
    repair = (package / "ensureWindowsTriangleRuntime.m").read_text(encoding="utf-8")
    assert "phase_stamps_beta.ensureWindowsTriangleRuntime(stampsRoot)" in runtime
    for name in ("libstdc++-6.dll", "libgcc_s_seh-1.dll", "libwinpthread-1.dll"):
        assert name in repair
    assert "copyfile(source,destination)" in repair
    assert '" -e phase-probe.1.node' in repair
    assert "phase-probe.2.edge" in repair


def test_train_utc_prefers_snap_master_over_diff_par(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    detector = (package / "autoDetectConfig.m").read_text(encoding="utf-8")
    schema = (package / "schema.m").read_text(encoding="utf-8")
    assert "masterTime = masterAcquisitionUtc(preprocFolder, cfg.master_date)" in detector
    assert "<MDATTR name=\"first_line_time\"" in detector
    assert "if ~isempty(masterTime)" in detector
    assert "cfg.utc_time = masterTime" in detector
    assert "verify it before using TRAIN" in detector
    assert "used by TRAIN/GACOS, not by StaMPS Step 8" in schema


def test_unix_installers_probe_native_triangle(phase_root):
    for relative in (
        "installer/prepare-macos-runtime.py",
        "installer/prepare-linux-runtime.py",
        "installer/install-phase-unix.py",
    ):
        source = (phase_root / relative).read_text(encoding="utf-8")
        assert '"-e"' in source
        assert "phase-probe.2.edge" in source
        assert "returncode != 0" in source


def test_stamps_step_8_can_be_retried_without_repeating_earlier_steps(phase_root):
    package = phase_root / "PHASE_Preprocessing" / "+phase_stamps_beta"
    schema = (package / "schema.m").read_text(encoding="utf-8")
    validation = (package / "validateConfig.m").read_text(encoding="utf-8")
    assert "{'1','2','3','4','5','6','7','8'}" in schema
    assert "{'1','2','3','4','5','6','7','8'}" in validation
