import shutil
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.requires_matlab
def test_matlab_can_open_repacked_mlapp(phase_root: Path, tmp_path: Path):
    """After tools/mlapp_roundtrip.py rewrites the file, MATLAB must still
    parse it (verifies canonical Office Open XML file ordering)."""
    src = phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp"
    dst = tmp_path / "PHASE_StaMPS.mlapp"
    shutil.copy(src, dst)
    # Trivial edit: identity transform (re-pack without content change)
    subprocess.check_call([
        sys.executable, str(phase_root / "tools/mlapp_roundtrip.py"),
        str(dst), r"(?s).*", r"\g<0>"  # match-everything, replace with same
    ])
    # Ask MATLAB to load the file
    cmd = [
        "matlab", "-batch",
        f"appdesigner.internal.serialization.app.readMLAPPFile('{dst}'); exit(0);"
    ]
    result = subprocess.run(cmd, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stderr.decode()
