#!/usr/bin/env python3
"""Build an Apple Silicon PHASE installer DMG from this checkout."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
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
    parser.add_argument("--gawk-source", type=Path,
                        help="Matching GNU awk 5.4.0 source archive for a public runtime DMG")
    parser.add_argument("--snaphu-source", type=Path,
                        help="Matching SNAPHU 2.0.7 source archive for a public runtime DMG")
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
    gawk_source = args.gawk_source.expanduser().resolve() if args.gawk_source else None
    snaphu_source = args.snaphu_source.expanduser().resolve() if args.snaphu_source else None
    build_tag = installer.release_tag(root)
    if runtime and build_tag and (not gawk_source or not snaphu_source):
        raise RuntimeError("Tagged runtime DMG requires GNU awk and SNAPHU source archives.")
    if gawk_source:
        if not runtime or gawk_source.name != "gawk-5.4.0.tar.xz":
            raise RuntimeError("GNU awk source requires the runtime and gawk-5.4.0.tar.xz.")
        digest = hashlib.sha256(gawk_source.read_bytes()).hexdigest()
        if digest != "3dd430f0cd3b4428c6c3f6afc021b9cd3c1f8c93f7a688dc268ca428a90b4ac1":
            raise RuntimeError("GNU awk source archive SHA-256 mismatch.")
        gawk_binary = runtime / "StaMPS" / "external" / "gawk" / "bin" / "gawk"
        result = subprocess.run([str(gawk_binary), "--version"],
                                capture_output=True, text=True, check=True)
        if not result.stdout.startswith("GNU Awk 5.4.0,"):
            raise RuntimeError("Bundled GNU awk binary does not match the source archive.")
    if snaphu_source:
        if not runtime or snaphu_source.name != "snaphu-v2.0.7.tar.gz":
            raise RuntimeError("SNAPHU source requires the runtime and snaphu-v2.0.7.tar.gz.")
        digest = hashlib.sha256(snaphu_source.read_bytes()).hexdigest()
        if digest != "c03ac126f9a964321bb5d6fb5b4004728368da268d7cb8407bb295a8abe5b262":
            raise RuntimeError("SNAPHU source archive SHA-256 mismatch.")
        snaphu_binary = runtime / "StaMPS" / "external" / "snaphu" / "bin" / "snaphu"
        result = subprocess.run([str(snaphu_binary), "-h"],
                                capture_output=True, text=True, check=False)
        if "snaphu v2.0.7" not in result.stdout + result.stderr:
            raise RuntimeError("Bundled SNAPHU binary does not match the source archive.")
    if args.sign_identity and runtime is None:
        raise RuntimeError("Signed PHASE release requires the complete StaMPS/TRAIN runtime.")
    if runtime:
        missing = installer.macos_psi_missing(runtime / "StaMPS", runtime / "TRAIN")
        if missing:
            raise RuntimeError("Incomplete Apple Silicon PSI runtime: " + ", ".join(missing))
        manifest = runtime / "phase-runtime.json"
        if args.sign_identity and not manifest.is_file():
            raise RuntimeError("Signed release requires a pinned phase-runtime.json manifest.")
        if manifest.is_file():
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            if (metadata.get("stampsCommit") != installer.STAMPS_COMMIT or
                    metadata.get("trainCommit") != installer.TRAIN_COMMIT or
                    metadata.get("platform") != "macos-arm64"):
                raise RuntimeError("macOS runtime revisions do not match PHASE release pins.")
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
        source = installer_dir / "PHASEInstaller.swift"
        if build_tag:
            generated_source = stage / "PHASEInstaller.swift"
            generated_source.write_text(
                source.read_text(encoding="utf-8").replace("v7.0.0 preview", build_tag),
                encoding="utf-8")
            source = generated_source
        build_env = os.environ.copy()
        if "DEVELOPER_DIR" not in build_env and Path("/Applications/Xcode.app/Contents/Developer").is_dir():
            build_env["DEVELOPER_DIR"] = "/Applications/Xcode.app/Contents/Developer"
        subprocess.run([
            "xcrun", "swiftc", "-parse-as-library", "-O", "-target", "arm64-apple-macos13.0",
            "-module-cache-path", str(stage / "swift-module-cache"),
            "-o", str(executable), str(source)
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
            "CFBundleVersion": build_tag[1:] if build_tag else "7.0.0",
            "CFBundleShortVersionString": build_tag[1:] if build_tag else "7.0.0",
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
        installer.embed_release_metadata(contents / "engine", root)
        if runtime:
            shutil.copytree(runtime / "StaMPS", contents / "StaMPS")
            shutil.copytree(runtime / "TRAIN", contents / "TRAIN")
        if gawk_source and snaphu_source:
            sources = volume / "Third-party sources"
            sources.mkdir()
            shutil.copy2(gawk_source, sources / gawk_source.name)
            shutil.copy2(snaphu_source, sources / snaphu_source.name)
            shutil.copy2(root / "Documenti" / "PHASE_7_THIRD_PARTY.md",
                         sources / "README.md")
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
