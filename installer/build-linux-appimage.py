#!/usr/bin/env python3
"""Build a PHASE x86_64 or aarch64 GUI installer AppImage on Linux."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "aarch64"}:
        raise RuntimeError("Build the AppImage on a supported Linux architecture.")
    appimagetool = shutil.which("appimagetool")
    if not appimagetool:
        raise RuntimeError("Install appimagetool before building the AppImage.")
    root = Path(__file__).resolve().parent.parent
    installer_dir = root / "installer"
    spec = importlib.util.spec_from_file_location(
        "phase_installer", installer_dir / "install-phase-unix.py")
    installer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(installer)
    installer.validate_runtime(root)
    output = args.output.expanduser().resolve()
    if root in output.parents:
        raise RuntimeError("Build the AppImage outside the source checkout.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="phase-appimage-") as temporary:
        appdir = Path(temporary) / "PHASE-Installer.AppDir"
        appdir.mkdir()
        launcher = appdir / "AppRun"
        launcher.write_text(
            '#!/bin/sh\nset -eu\nappdir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)\n'
            'exec "$appdir/usr/bin/phase-installer" "$@"\n', encoding="utf-8"
        )
        launcher.chmod(0o755)
        (appdir / "phase-installer.desktop").write_text(
            "[Desktop Entry]\nType=Application\nName=PHASE Installer\n"
            "Exec=phase-installer\nIcon=phase-installer\nTerminal=false\n"
            "Categories=Science;Education;\n", encoding="utf-8"
        )
        shutil.copy2(root / "Logo_square.png", appdir / "phase-installer.png")
        (appdir / ".DirIcon").symlink_to("phase-installer.png")
        bin_dir = appdir / "usr" / "bin"
        bin_dir.mkdir(parents=True)
        shutil.copy2(installer_dir / "PHASE-Linux-Installer.sh", bin_dir / "phase-installer")
        (bin_dir / "phase-installer").chmod(0o755)
        resources = appdir / "usr" / "share" / "phase"
        resources.mkdir(parents=True)
        shutil.copy2(installer_dir / "install-phase-unix.py", resources)
        installer.SOURCE_FOR_COPY = root
        shutil.copytree(root, resources / "engine", ignore=installer.ignored)
        env = dict(os.environ, ARCH=platform.machine())
        subprocess.run([appimagetool, str(appdir), str(output)], check=True, env=env)
    print(output)


if __name__ == "__main__":
    main()
