"""Verify PHASE's setupPythonEnvironment writes %APPDATA%\\PHASE\\python.txt
so StaMPS's Windows .bat shim can honor the same interpreter.

This test invokes MATLAB, so it is marked requires_matlab + windows_only and
will skip on Linux/macOS or where MATLAB is unavailable.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


@pytest.mark.windows_only
@pytest.mark.requires_matlab
def test_phase_writes_python_txt(phase_root: Path, tmp_path: Path,
                                 monkeypatch: pytest.MonkeyPatch) -> None:
    """Run setupPythonEnvironment and verify python.txt is created."""
    monkeypatch.setenv("APPDATA", str(tmp_path))
    cmd = [
        "matlab", "-batch",
        f"cd('{phase_root}/MatlabFunctions'); setupPythonEnvironment; exit(0);",
    ]
    proc = subprocess.run(cmd, capture_output=True, timeout=120)
    assert proc.returncode == 0, proc.stderr.decode(errors="replace")
    pt = tmp_path / "PHASE" / "python.txt"
    assert pt.exists(), "python.txt not created by PHASE"
    content = pt.read_text(encoding="utf-8").strip()
    assert content.endswith("python.exe") or content.endswith("python3"), (
        f"unexpected interpreter path: {content!r}"
    )
