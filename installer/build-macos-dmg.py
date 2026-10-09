#!/usr/bin/env python3
"""Build an Apple Silicon PHASE installer DMG from this checkout."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
import platform
import plistlib
import shutil
import subprocess
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sign-identity", help="Developer ID Application identity")
    parser.add_argument("--notary-profile", help="notarytool keychain profile")
    parser.add_argument("--runtime", type=Path,
                        help="Prepared StaMPS/TRAIN runtime for a local PSI installer preview")
    args = parser.parse_args()
    if args.notary_profile and not args.sign_identity:
        parser.error("--notary-profile requires --sign-identity")
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise RuntimeError("Build the macOS installer on an Apple Silicon Mac.")
    root = Path(__file__).resolve().parent.parent
    installer_dir = root / "installer"
    source_spec = importlib.util.spec_from_file_location(
        "phase_installer", installer_dir / "install-phase-unix.py")
    installer = importlib.util.module_from_spec(source_spec)
    source_spec.loader.exec_module(installer)
    installer.validate_runtime(root)
    runtime = args.runtime.expanduser().resolve() if args.runtime else None
    if runtime:
        missing = installer.macos_psi_missing(runtime / "StaMPS", runtime / "TRAIN")
        if missing:
            raise RuntimeError("Incomplete Apple Silicon PSI runtime: " + ", ".join(missing))
    output = args.output.expanduser().resolve()
    if root in output.parents:
        raise RuntimeError("Build the DMG outside the source checkout.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="phase-macos-dmg-") as temporary:
        stage = Path(temporary)
        volume = stage / "volume"
        volume.mkdir()
        app = volume / "PHASE Installer.app"
        bundle = app / "Contents"
        contents = bundle / "Resources"
        executable = bundle / "MacOS" / "PHASEInstaller"
        contents.mkdir(parents=True)
        executable.parent.mkdir(parents=True)
        build_env = os.environ.copy()
        if "DEVELOPER_DIR" not in build_env and Path("/Applications/Xcode.app/Contents/Developer").is_dir():
            build_env["DEVELOPER_DIR"] = "/Applications/Xcode.app/Contents/Developer"
        subprocess.run([
            "xcrun", "swiftc", "-parse-as-library", "-O", "-target", "arm64-apple-macos13.0",
            "-module-cache-path", str(stage / "swift-module-cache"),
            "-o", str(executable), str(installer_dir / "PHASEInstaller.swift")
        ], env=build_env, check=True)
        shutil.copy2(installer_dir / "install-phase-unix.py", contents)
        shutil.copy2(installer_dir / "find-python-macos.sh", contents)
        shutil.copy2(installer_dir / "PHASE.icns", contents)
        shutil.copy2(root / "PHASE_logo.png", contents)
        plist_path = bundle / "Info.plist"
        plist = {
            "CFBundleName": "PHASE Installer",
            "CFBundleDisplayName": "PHASE Installer",
            "CFBundleIdentifier": "org.phaseinsar.phase.installer",
            "CFBundleVersion": "7.0.0",
            "CFBundleShortVersionString": "7.0.0",
            "CFBundleExecutable": "PHASEInstaller",
            "CFBundlePackageType": "APPL",
            "CFBundleIconFile": "PHASE.icns",
            "LSMinimumSystemVersion": "13.0",
            "NSHighResolutionCapable": True,
        }
        with plist_path.open("wb") as stream:
            plistlib.dump(plist, stream)
        installer.SOURCE_FOR_COPY = root
        shutil.copytree(root, contents / "engine", ignore=installer.ignored)
        if runtime:
            shutil.copytree(runtime / "StaMPS", contents / "StaMPS")
            shutil.copytree(runtime / "TRAIN", contents / "TRAIN")
        signing = ["codesign", "--force", "--deep", "--sign",
                   args.sign_identity or "-"]
        if args.sign_identity:
            signing.extend(["--options", "runtime", "--timestamp"])
        subprocess.run([*signing, str(app)], check=True)
        subprocess.run([
            "hdiutil", "create", "-volname", "PHASE 7 Installer",
            "-srcfolder", str(volume), "-format", "UDZO", "-ov", str(output)
        ], check=True)
        if args.notary_profile:
            subprocess.run([
                "xcrun", "notarytool", "submit", str(output),
                "--keychain-profile", args.notary_profile, "--wait"
            ], check=True)
            subprocess.run(["xcrun", "stapler", "staple", str(output)], check=True)
    print(output)


if __name__ == "__main__":
    main()
