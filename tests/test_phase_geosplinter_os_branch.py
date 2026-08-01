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
