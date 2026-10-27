"""Exercise the installed updater as a separate process, including on Windows."""

import json
from pathlib import Path
import shutil
import subprocess
import sys

import phase_update
import pytest


def _engine(folder: Path, tag: str) -> None:
    folder.mkdir(parents=True)
    for relative in phase_update.REQUIRED_FILES:
        target = folder / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("fixture", encoding="utf-8")
    (folder / "phase-release.json").write_text(
        json.dumps({"tag": tag, "updateSchema": 1}), encoding="utf-8")


@pytest.mark.parametrize("partial_runtime", [False, True])
def test_apply_replaces_engine_when_updater_runs_from_it(tmp_path, partial_runtime):
    prefix = tmp_path / "PHASE with spaces"
    engine = prefix / "engine"
    _engine(engine, "v7.0.3")
    shutil.copy2(Path(phase_update.__file__), engine / "phase_update.py")
    (engine / "StaMPS").mkdir()
    (engine / "StaMPS" / "runtime.txt").write_text("keep", encoding="utf-8")
    (engine / "TRAIN").mkdir()
    (engine / "TRAIN" / "runtime.txt").write_text("keep", encoding="utf-8")
    (prefix / "install.json").write_text(
        json.dumps({"updateSchema": 1, "version": "v7.0.3"}), encoding="utf-8")
    pending = prefix / "pending-update"
    _engine(pending / "engine", "v7.0.4")
    if partial_runtime:
        for name in ("StaMPS", "TRAIN"):
            staged = pending / "engine" / name
            staged.mkdir()
            (staged / "partial.txt").write_text("incomplete", encoding="utf-8")
    (pending / "update.json").write_text(
        json.dumps({"updateSchema": 1, "version": "v7.0.4"}), encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(engine / "phase_update.py"), "apply", "--prefix", str(prefix)],
        text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert phase_update.read_json(prefix / "install.json")["version"] == "v7.0.4"
    assert (prefix / "engine" / "StaMPS" / "runtime.txt").read_text() == "keep"
    assert (prefix / "engine" / "TRAIN" / "runtime.txt").read_text() == "keep"
    assert not (prefix / "engine" / "StaMPS" / "partial.txt").exists()
    assert not (prefix / "engine" / "TRAIN" / "partial.txt").exists()
