#!/usr/bin/env python3
"""Build the phase7-engine.zip asset for a PHASE 7 GitHub release."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True, help="Release tag, e.g. v7.0.0")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--dev-build", action="store_true",
                        help="Skip clean tagged checkout requirement for local trials")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    updater_spec = importlib.util.spec_from_file_location("phase_updater", root / "phase_update.py")
    updater = importlib.util.module_from_spec(updater_spec)
    updater_spec.loader.exec_module(updater)
    updater.version_tuple(args.tag)
    if not args.dev_build:
        exact_tag = subprocess.check_output(
            ["git", "describe", "--tags", "--exact-match"], cwd=root, text=True
        ).strip()
        if exact_tag != args.tag:
            raise RuntimeError(f"Checkout tag {exact_tag} does not match {args.tag}")
        if subprocess.check_output(["git", "status", "--porcelain"], cwd=root):
            raise RuntimeError("Build release assets from a clean checkout.")
    installer_spec = importlib.util.spec_from_file_location(
        "phase_installer", root / "installer" / "install-phase-unix.py")
    installer = importlib.util.module_from_spec(installer_spec)
    installer_spec.loader.exec_module(installer)
    installer.validate_runtime(root)
    output = args.output.expanduser().resolve()
    if output.name != updater.ASSET_NAME:
        raise RuntimeError(f"Release asset must be named {updater.ASSET_NAME}")
    if root in output.parents:
        raise RuntimeError("Build the release asset outside the source checkout.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="phase-release-") as temporary:
        engine = Path(temporary) / "engine"
        installer.SOURCE_FOR_COPY = root
        shutil.copytree(root, engine, ignore=installer.ignored)
        (engine / "phase-release.json").write_text(json.dumps({
            "tag": args.tag, "updateSchema": 1
        }, indent=2) + "\n", encoding="utf-8")
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for file in sorted(engine.rglob("*")):
                if file.is_file():
                    archive.write(file, file.relative_to(engine).as_posix())
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    print(json.dumps({"asset": str(output), "sha256": digest}, indent=2))


if __name__ == "__main__":
    main()
