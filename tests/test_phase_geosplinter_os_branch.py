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
def test_geosplinter_has_isunix_branch(phase_root: Path, filename: str):
    path = phase_root / "MatlabFunctions" / filename
    text = path.read_text(encoding="utf-8")
    assert "if isunix" in text, f"{filename}: no isunix branch found"
    assert "tempname()" in text, f"{filename}: missing Windows tempname branch"
    # Ensure no raw `sprintf('%s < %s', ...)` remains without an isunix guard
    # (the raw pattern may legitimately appear INSIDE the isunix branch).
    lines = text.splitlines()
    in_isunix = False
    for i, line in enumerate(lines):
        if "if isunix" in line: in_isunix = True
        elif line.strip() == "end" and in_isunix: in_isunix = False
        elif "sprintf(" in line and "'%s < %s'" in line and not in_isunix:
            pytest.fail(f"{filename}:{i+1}: raw stdin redirect outside isunix guard")
