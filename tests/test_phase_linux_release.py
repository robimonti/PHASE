"""Linux installer release guards (host-independent checks)."""

import importlib.util
import subprocess


def _load(root):
    path = root / "installer" / "prepare-linux-runtime.py"
    spec = importlib.util.spec_from_file_location("phase_linux_runtime", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_linux_runtime_uses_tested_fork_revisions(phase_root):
    module = _load(phase_root)
    assert module.STAMPS_REPO == "https://github.com/pyccino/StaMPS.git"
    assert module.TRAIN_REPO == "https://github.com/pyccino/TRAIN.git"
    assert module.STAMPS_COMMIT == "7cabf05eddf8ebe8694e5346fe0f9d48aaef4962"
    assert module.TRAIN_COMMIT == "6d0273ae67d2a9f07a696b6a14298ef2c31607d8"
    windows = (phase_root / "installer" / "install-phase.ps1").read_text(encoding="utf-8-sig")
    assert module.STAMPS_COMMIT in windows
    assert module.TRAIN_COMMIT in windows
    assert "-Commit $Script:StampsCommit" in windows
    assert "-Commit $Script:TrainCommit" in windows


def test_linux_appimage_contains_native_runtime_preparer(phase_root):
    builder = (phase_root / "installer" / "build-linux-appimage.py").read_text()
    wizard = (phase_root / "installer" / "PHASE-Linux-Installer.sh").read_text()
    assert "prepare-linux-runtime.py" in builder
    assert "prepare-linux-runtime.py" in wizard
    assert '--stamps "$runtime_tmp/runtime/StaMPS"' in wizard
    assert '--train "$runtime_tmp/runtime/TRAIN"' in wizard
    subprocess.run(["sh", "-n", str(phase_root / "installer" / "PHASE-Linux-Installer.sh")],
                   check=True)
