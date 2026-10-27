import pytest
from pathlib import Path


GEOSPLINTER_FILES = [
    "STmodel_DET1D.m",
    "STmodel_DET2D.m",
    "STmodel_STC1D.m",
    "STmodel_STC2D.m",
    "ModellingInTime.m",
]


@pytest.mark.parametrize("filename", GEOSPLINTER_FILES)
def test_geosplinter_uses_cross_platform_direct_runner(phase_root: Path, filename: str):
    path = phase_root / "MatlabFunctions" / filename
    text = path.read_text(encoding="utf-8")
    assert "runGeoSplinter(" in text, f"{filename}: direct runner not used"
    assert "tempname()" not in text, f"{filename}: legacy Windows batch file remains"
    assert "'%s < %s'" not in text, f"{filename}: shell stdin redirection remains"
    assert "system(job_execution" not in text


def test_direct_runner_resolves_windows_executable_and_redirects_job(phase_root: Path):
    runner = (phase_root / "MatlabFunctions" / "runGeoSplinter.m").read_text(
        encoding="utf-8"
    )
    assert "java.lang.ProcessBuilder" in runner
    assert "redirectInput(java.io.File(jobFile))" in runner
    assert "endsWith(lower(executable),'.exe')" in runner
    assert "candidate.getCanonicalPath()" in runner
    assert "workingDir = fileparts(fileparts(jobFile))" in runner


def test_model_jobs_support_absolute_project_paths(phase_root: Path):
    functions = phase_root / "MatlabFunctions"
    helper = (functions / "geoSplinterJobPath.m").read_text(encoding="utf-8")
    assert "folderParent" in helper and "jobParent" in helper
    assert "value = fullfile('.',folderName" in helper
    for name in ("jobFile_analysis.m", "jobFile_analysis2D.m", "jobFile_synthesis.m"):
        source = (functions / name).read_text(encoding="utf-8")
        assert "job_file = fullfile(gS_job_path" in source
        assert "geoSplinterJobPath(" in source
        assert "fullfile('.', gS_job_path" not in source
    for name in GEOSPLINTER_FILES:
        source = (functions / name).read_text(encoding="utf-8")
        assert "fullfile('.', gS_" not in source
    smoke = (phase_root / "+phase_model_beta" / "geoSplinterSelfTest.m").read_text(
        encoding="utf-8"
    )
    assert "absoluteRoot = [tempname ' PHASE model test']" in smoke
