#!/usr/bin/env python3
"""Resumable ASF downloader with atomic progress reporting for PHASE."""

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path


CHUNK_SIZE = 1024 * 1024
PROGRESS_INTERVAL = 0.25
MAX_ATTEMPTS = 3


class StopRequested(Exception):
    pass


def atomic_json(path, payload):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, destination)


def load_json(path):
    with open(path, "r", encoding="utf-8") as stream:
        return json.load(stream)


def safe_name(raw):
    name = os.path.basename(str(raw or "").strip())
    if not name:
        raise ValueError("A download entry has no file name.")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise ValueError(f"Unsafe download file name: {name}")
    if not name.lower().endswith(".zip"):
        name += ".zip"
    return name


def normalize_entries(raw):
    if not isinstance(raw, list) or not raw:
        raise ValueError("The download manifest contains no files.")
    entries = []
    seen = set()
    for item in raw:
        if not isinstance(item, dict):
            raise ValueError("Every download manifest entry must be an object.")
        name = safe_name(item.get("name") or item.get("sceneName"))
        url = str(item.get("url") or "").strip()
        if not url.lower().startswith("https://"):
            raise ValueError(f"Invalid HTTPS download URL for {name}.")
        if name in seen:
            continue
        seen.add(name)
        try:
            expected = max(0, int(float(item.get("sizeBytes") or 0)))
        except (TypeError, ValueError):
            expected = 0
        entries.append({"name": name, "url": url, "sizeBytes": expected})
    if not entries:
        raise ValueError("The download manifest contains no unique files.")
    return entries


def authenticated_session(credentials_path):
    try:
        import asf_search as asf
    except ImportError as exc:
        raise RuntimeError(
            "Python package 'asf_search' is missing. Install it in the Python "
            "environment selected by PHASE."
        ) from exc
    credentials = load_json(credentials_path)
    username = str(credentials.get("username") or "").strip()
    password = str(credentials.get("password") or "")
    if not username or not password:
        raise RuntimeError("Saved Earthdata credentials are incomplete.")
    return asf.ASFSession().auth_with_creds(username, password)


class Progress:
    def __init__(self, path, entries, kind):
        self.path = Path(path)
        self.entries = entries
        self.kind = kind
        self.current_index = 0
        self.current_file = ""
        self.current_bytes = 0
        self.current_total = 0
        self.completed_bytes = 0
        self.completed_files = 0
        self.skipped_files = 0
        self.failed_files = 0
        self.message = "Preparing download…"
        self.phase = "preparing"
        self.last_write = 0
        self.expected_total = sum(item["sizeBytes"] for item in entries)

    def percentage(self):
        if self.expected_total > 0 and all(item["sizeBytes"] > 0 for item in self.entries):
            value = 100 * (self.completed_bytes + self.current_bytes) / self.expected_total
        else:
            fraction = 0
            if self.current_total > 0:
                fraction = min(1, self.current_bytes / self.current_total)
            value = 100 * (self.completed_files + fraction) / max(1, len(self.entries))
        return round(max(0, min(100, value)), 2)

    def payload(self):
        return {
            "kind": self.kind,
            "phase": self.phase,
            "active": self.phase in {"preparing", "downloading", "retrying", "stopping"},
            "percentage": self.percentage(),
            "currentIndex": self.current_index,
            "totalFiles": len(self.entries),
            "currentFile": self.current_file,
            "currentBytes": self.current_bytes,
            "currentTotalBytes": self.current_total,
            "completedBytes": self.completed_bytes,
            "expectedTotalBytes": self.expected_total,
            "completedFiles": self.completed_files,
            "skippedFiles": self.skipped_files,
            "failedFiles": self.failed_files,
            "message": self.message,
            "pid": os.getpid(),
            "updatedAt": time.time(),
        }

    def write(self, force=False):
        now = time.monotonic()
        if force or now - self.last_write >= PROGRESS_INTERVAL:
            atomic_json(self.path, self.payload())
            self.last_write = now


def stop_if_requested(stop_path):
    if stop_path and Path(stop_path).exists():
        raise StopRequested("Download stopped by the user.")


def response_total(response, offset):
    content_range = response.headers.get("Content-Range", "")
    match = re.search(r"/(\d+)$", content_range)
    if match:
        return int(match.group(1))
    length = response.headers.get("Content-Length")
    if length:
        return offset + int(length)
    return 0


def download_one(session, entry, destination, stop_path, progress):
    target = destination / entry["name"]
    partial = target.with_name(target.name + ".part")
    expected = entry["sizeBytes"]

    if target.is_file():
        actual = target.stat().st_size
        if expected == 0 or actual == expected:
            return "skipped", actual
        target.unlink()
    if partial.is_file() and expected:
        partial_size = partial.stat().st_size
        if partial_size == expected:
            os.replace(partial, target)
            return "downloaded", partial_size
        if partial_size > expected:
            partial.unlink()

    for attempt in range(1, MAX_ATTEMPTS + 1):
        stop_if_requested(stop_path)
        offset = partial.stat().st_size if partial.is_file() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            with session.get(
                entry["url"], stream=True, headers=headers, timeout=(30, 90)
            ) as response:
                response.raise_for_status()
                append = offset > 0 and response.status_code == 206
                if offset and not append:
                    offset = 0
                total = response_total(response, offset) or expected
                progress.current_bytes = offset
                progress.current_total = total
                progress.phase = "downloading"
                progress.message = (
                    f"Downloading image {progress.current_index} of "
                    f"{len(progress.entries)}"
                )
                progress.write(force=True)

                mode = "ab" if append else "wb"
                with partial.open(mode) as stream:
                    for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                        stop_if_requested(stop_path)
                        if not chunk:
                            continue
                        stream.write(chunk)
                        progress.current_bytes += len(chunk)
                        progress.write()
                    stream.flush()
                    os.fsync(stream.fileno())

            actual = partial.stat().st_size
            verified = total or expected
            if verified and actual != verified:
                raise IOError(
                    f"Incomplete file {entry['name']}: {actual} of {verified} bytes."
                )
            os.replace(partial, target)
            return "downloaded", actual
        except StopRequested:
            raise
        except Exception as exc:
            if attempt >= MAX_ATTEMPTS:
                raise RuntimeError(
                    f"{entry['name']} failed after {MAX_ATTEMPTS} attempts: {exc}"
                ) from exc
            progress.phase = "retrying"
            progress.message = (
                f"Retrying image {progress.current_index} of {len(progress.entries)} "
                f"(attempt {attempt + 1}/{MAX_ATTEMPTS})"
            )
            progress.write(force=True)
            time.sleep(min(2**attempt, 5))
    raise RuntimeError(f"Download failed for {entry['name']}.")


def run(args, session=None):
    manifest = load_json(args.manifest)
    entries = normalize_entries(manifest.get("files", manifest))
    destination = Path(args.destination)
    destination.mkdir(parents=True, exist_ok=True)
    progress = Progress(args.progress, entries, args.kind)
    result = {
        "kind": args.kind,
        "status": "running",
        "downloaded": [],
        "skipped": [],
        "failed": [],
    }
    progress.write(force=True)
    if session is None:
        session = authenticated_session(args.credentials)

    try:
        for index, entry in enumerate(entries, 1):
            stop_if_requested(args.stop)
            progress.current_index = index
            progress.current_file = entry["name"]
            progress.current_bytes = 0
            progress.current_total = entry["sizeBytes"]
            progress.write(force=True)
            try:
                status, size = download_one(
                    session, entry, destination, args.stop, progress
                )
                if status == "skipped":
                    progress.skipped_files += 1
                    result["skipped"].append(entry["name"])
                else:
                    result["downloaded"].append(entry["name"])
                progress.completed_files += 1
                progress.completed_bytes += size
                progress.current_bytes = 0
                progress.current_total = 0
                progress.message = (
                    f"Completed {progress.completed_files} of {len(entries)} images"
                )
                progress.write(force=True)
            except StopRequested:
                raise
            except Exception as exc:
                progress.failed_files += 1
                result["failed"].append(
                    {"name": entry["name"], "message": str(exc)}
                )
                progress.current_bytes = 0
                progress.current_total = 0
                progress.message = str(exc)
                progress.write(force=True)

        result["status"] = "failed" if result["failed"] else "completed"
        progress.phase = result["status"]
        progress.message = (
            f"Completed {progress.completed_files} of {len(entries)} images"
            if not result["failed"]
            else f"{len(result['failed'])} image download(s) failed"
        )
    except StopRequested as exc:
        result["status"] = "stopped"
        progress.phase = "stopped"
        progress.message = str(exc)

    result.update(
        {
            "totalFiles": len(entries),
            "completedFiles": progress.completed_files,
            "downloadedCount": len(result["downloaded"]),
            "skippedCount": len(result["skipped"]),
            "failedCount": len(result["failed"]),
        }
    )
    progress.write(force=True)
    atomic_json(args.output, result)
    return result


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--credentials", required=True)
    parser.add_argument("--progress", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stop", required=True)
    parser.add_argument("--kind", choices=("initial", "update"), required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        result = run(args)
    except Exception as exc:
        atomic_json(
            args.output,
            {
                "kind": args.kind,
                "status": "failed",
                "downloaded": [],
                "skipped": [],
                "failed": [{"name": "", "message": str(exc)}],
                "downloadedCount": 0,
                "skippedCount": 0,
                "failedCount": 1,
            },
        )
        print(str(exc), file=sys.stderr)
        return 1
    if result["status"] == "failed":
        return 1
    if result["status"] == "stopped":
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
