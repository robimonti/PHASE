import subprocess
from pathlib import Path

import pytest


@pytest.mark.requires_matlab
@pytest.mark.nightly
def test_phase_preprocessing_launches(phase_root: Path):
    """Minimal assertion: PHASE_Preprocessing.mlapp can be loaded without error."""
    cmd = [
        "matlab", "-batch",
        f"cd('{phase_root}'); "
        f"appdesigner.internal.serialization.app.readMLAPPFile("
        f"'PHASE_Preprocessing.mlapp'); exit(0);"
    ]
    result = subprocess.run(cmd, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stderr.decode()
