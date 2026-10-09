#!/usr/bin/env python3
"""Stage PHASE 7 GitHub releases and apply them before MATLAB starts.

Only the engine is replaced. Projects live outside the installation; StaMPS
and TRAIN bundled by the installer are carried forward into the new engine.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import sys
import tempfile
import urllib.request
import zipfile


API_URL = "https://api.github.com/repos/robimonti/PHASE/releases/latest"
ASSET_NAME = "phase7-engine.zip"
REQUIRED_FILES = (
    "PHASE_Hub.m", "PHASE_Hub_UI.html", "PHASE_logo.png",
    "PHASE_mod1a.png", "PHASE_mod1b.png", "PHASE_mod2.png",
    "+phase_hub/App.m", "+phase_project/open.m",
    "+phase_project/workflowStatus.m",
    "PHASE_Preprocessing_beta.m", "+phase_model_beta/App.m",
    "PHASE_Preprocessing/+phase_stamps_beta/App.m",
)
RELEASE_TAG = re.compile(r"^v(\d+)\.(\d+)\.(\d+)$")
MAX_ARCHIVE_BYTES = 1_000_000_000


def read_json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise RuntimeError(f"Invalid JSON object: {path}")
    return data


def managed_prefix(path: str) -> tuple[Path, dict]:
    prefix = Path(path).expanduser().resolve()
    metadata = prefix / "install.json"
    if not metadata.is_file() or not (prefix / "engine" / "PHASE_Hub.m").is_file():
        raise RuntimeError("PHASE update requires a managed PHASE 7 installation.")
    info = read_json(metadata)
    if info.get("updateSchema") != 1:
        raise RuntimeError("This installation does not support in-app updates; reinstall PHASE 7.")
    return prefix, info


def version_tuple(tag: str) -> tuple[int, int, int]:
    match = RELEASE_TAG.fullmatch(tag)
    if not match or int(match.group(1)) < 7:
        raise RuntimeError(f"Unsupported PHASE release tag: {tag}")
    return tuple(int(part) for part in match.groups())


def release_info() -> dict:
    request = urllib.request.Request(API_URL, headers={
        "Accept": "application/vnd.github+json", "User-Agent": "PHASE-7-Updater"
    })
    with urllib.request.urlopen(request, timeout=15) as response:
        data = json.load(response)
    if data.get("draft") or data.get("prerelease"):
        raise RuntimeError("Latest GitHub release is not a stable PHASE release.")
    tag = data.get("tag_name", "")
    version_tuple(tag)
    asset = next((item for item in data.get("assets", [])
                  if item.get("name") == ASSET_NAME), None)
    if asset is None:
        raise RuntimeError(f"Release {tag} has no {ASSET_NAME} asset.")
    digest = asset.get("digest", "")
    if not re.fullmatch(r"sha256:[0-9a-fA-F]{64}", digest):
        raise RuntimeError(f"Release {tag} has no SHA-256 asset digest.")
    url = asset.get("browser_download_url", "")
    expected = f"https://github.com/robimonti/PHASE/releases/download/{tag}/{ASSET_NAME}"
    if url != expected:
        raise RuntimeError("Unexpected PHASE release asset URL.")
    return {"version": tag, "url": url, "digest": digest.lower()}


def check(prefix: Path, installed: dict) -> dict:
    current = installed.get("version", "dev")
    try:
        available = release_info()
    except RuntimeError as error:
        if (str(error).startswith("Unsupported PHASE release tag:") or
                "has no phase7-engine.zip asset" in str(error)):
            return {"current": current, "available": None, "updateAvailable": False}
        raise
    newer = current == "dev" or version_tuple(available["version"]) > version_tuple(current)
    return {"current": current, "available": available["version"], "updateAvailable": newer}


def download_verified(url: str, digest: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "PHASE-7-Updater"})
    hash_value = hashlib.sha256()
    total = 0
    with urllib.request.urlopen(request, timeout=60) as response, destination.open("wb") as output:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > MAX_ARCHIVE_BYTES:
                raise RuntimeError("PHASE release archive is unexpectedly large.")
            hash_value.update(chunk)
            output.write(chunk)
    if hash_value.hexdigest() != digest.split(":", 1)[1]:
        raise RuntimeError("PHASE release checksum does not match GitHub metadata.")


def extract_engine(archive: Path, destination: Path, tag: str) -> None:
    destination.mkdir()
    total = 0
    with zipfile.ZipFile(archive) as package:
        for item in package.infolist():
            member = PurePosixPath(item.filename)
            if (member.is_absolute() or not member.parts or
                    any(part in {"", ".", ".."} for part in member.parts) or
                    "\\" in item.filename or ":" in member.parts[0]):
                raise RuntimeError(f"Unsafe path in PHASE release: {item.filename}")
            mode = (item.external_attr >> 16) & 0o170000
            if mode not in {0, stat.S_IFREG, stat.S_IFDIR}:
                raise RuntimeError(f"Unsupported archive entry: {item.filename}")
            total += item.file_size
            if total > MAX_ARCHIVE_BYTES:
                raise RuntimeError("PHASE release contents are unexpectedly large.")
            target = destination.joinpath(*member.parts)
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists():
                    raise RuntimeError(f"Duplicate PHASE release entry: {item.filename}")
                with package.open(item) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output)
    missing = [name for name in REQUIRED_FILES if not (destination / name).is_file()]
    if missing:
        raise RuntimeError("Incomplete PHASE 7 update: " + ", ".join(missing))
    manifest = read_json(destination / "phase-release.json")
    if manifest.get("tag") != tag or manifest.get("updateSchema") != 1:
        raise RuntimeError("PHASE release manifest does not match the GitHub tag.")


def prepare(prefix: Path, installed: dict) -> dict:
    release = release_info()
    current = installed.get("version", "dev")
    if current != "dev" and version_tuple(release["version"]) <= version_tuple(current):
        return {"prepared": False, "version": current, "reason": "already current"}
    with tempfile.TemporaryDirectory(prefix=".phase-download-", dir=prefix) as temporary:
        temporary = Path(temporary)
        archive = temporary / ASSET_NAME
        download_verified(release["url"], release["digest"], archive)
        engine = temporary / "engine"
        extract_engine(archive, engine, release["version"])
        staged = temporary / "ready"
        staged.mkdir()
        engine.rename(staged / "engine")
        (staged / "update.json").write_text(json.dumps({
            "updateSchema": 1, "version": release["version"],
            "digest": release["digest"],
        }, indent=2) + "\n", encoding="utf-8")
        pending = prefix / "pending-update"
        previous = temporary / "previous-pending"
        if pending.exists():
            if not (pending / "update.json").is_file():
                raise RuntimeError("Unrecognized pending-update directory; inspect it manually.")
            pending.rename(previous)
        try:
            staged.rename(pending)
        except Exception:
            if previous.exists():
                previous.rename(pending)
            raise
    return {"prepared": True, "version": release["version"],
            "message": "Close MATLAB, then launch PHASE again to apply the update."}


def apply(prefix: Path, installed: dict) -> dict:
    pending = prefix / "pending-update"
    if not pending.exists():
        return {"applied": False, "reason": "no pending update"}
    if not (pending / "update.json").is_file():
        raise RuntimeError("Unrecognized pending-update directory; inspect it manually.")
    update = read_json(pending / "update.json")
    tag = update.get("version", "")
    version_tuple(tag)
    candidate = pending / "engine"
    missing = [name for name in REQUIRED_FILES if not (candidate / name).is_file()]
    if missing or read_json(candidate / "phase-release.json").get("tag") != tag:
        raise RuntimeError("Pending PHASE update is incomplete.")
    old_engine = prefix / "engine"
    for name in ("StaMPS", "TRAIN"):
        source = old_engine / name
        if source.exists():
            target = candidate / name
            if target.exists():
                raise RuntimeError(f"PHASE release unexpectedly contains {name}.")
            shutil.copytree(source, target)
    for name in ("project.conf.template", "PHASE_Preprocessing/input_StaMPS.mat",
                 "PHASE_Preprocessing/input_preprocessing.mat"):
        source = old_engine / name
        target = candidate / name
        if source.is_file() and not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    backup = prefix / "backups" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup.mkdir(parents=True)
    moved_old = False
    moved_new = False
    try:
        old_engine.rename(backup / "engine")
        moved_old = True
        candidate.rename(old_engine)
        moved_new = True
        revised = dict(installed, version=tag, updatedAt=datetime.now(timezone.utc).isoformat())
        temporary_info = prefix / "install.json.new"
        temporary_info.write_text(json.dumps(revised, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary_info, prefix / "install.json")
    except Exception:
        if moved_new:
            old_engine.rename(candidate)
        if moved_old:
            (backup / "engine").rename(old_engine)
        raise
    shutil.rmtree(pending)
    return {"applied": True, "version": tag, "backup": str(backup / "engine")}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "prepare", "apply"))
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    try:
        prefix, installed = managed_prefix(args.prefix)
        result = {"check": check, "prepare": prepare, "apply": apply}[args.action](
            prefix, installed)
        print(json.dumps(result, ensure_ascii=False))
        return 0
    except Exception as error:
        print(json.dumps({"error": str(error)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
