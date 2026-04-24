"""Runtime test: executes the shipping guard block and verifies degradation."""
from __future__ import annotations
import re
import subprocess
import zipfile
from pathlib import Path

import pytest


def _extract_guard(mlapp_path: Path) -> str:
    with zipfile.ZipFile(mlapp_path) as z:
        xml = z.read("matlab/document.xml").decode("utf-8")
    # Anchored on the leading/trailing sentinels — robust against nested blocks.
    m = re.search(
        r"(% begin TRAIN availability check[\s\S]*?% end TRAIN availability check)",
        xml)
    assert m, "Guard block not found — Change #1 may not be applied"
    return m.group(1)


@pytest.mark.windows_only
@pytest.mark.requires_matlab
def test_degradation_when_train_missing(phase_root: Path, tmp_path: Path):
    guard = _extract_guard(
        phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")

    driver = tmp_path / "driver.m"
    driver.write_text(
        "warning('on', 'all');\n"       # defensive: ensure warnings are not suppressed
        "train_flag = 0;\n"
        "lastwarn('');\n"               # reset IMMEDIATELY before the guard so
                                        # only the guard's warning is captured
        + guard + "\n"
        "[msg, id] = lastwarn();\n"
        "fid = fopen('result.txt', 'w');\n"
        "fprintf(fid, 'train_flag=%d\\nid=%s\\n', train_flag, id);\n"
        "fclose(fid);\n",
        encoding="utf-8")

    # MATLAB's cd() on Windows accepts forward slashes; use as_posix() to
    # avoid any ambiguity with backslashes or single-quote escaping.
    mat_cwd = tmp_path.as_posix()
    result = subprocess.run(
        ["matlab", "-batch", f"cd('{mat_cwd}'); driver"],
        capture_output=True, timeout=120)
    assert result.returncode == 0, result.stderr.decode(errors="replace")

    contents = (tmp_path / "result.txt").read_text(encoding="utf-8")
    assert "train_flag=1" in contents, \
        f"Expected local train_flag downgrade to 1, got:\n{contents}"
    assert "id=StaMPS:phase:trainNotAvailable" in contents, \
        f"Expected canonical warning id in lastwarn, got:\n{contents}"
