#!/usr/bin/env python3
"""Build the tested StaMPS/TRAIN forks locally for a Linux PHASE install.

Third-party sources are fetched on the user's machine. SNAPHU, GNU awk and
the build toolchain remain system prerequisites; no third-party binary is
redistributed by this helper.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


STAMPS_REPO = "https://github.com/pyccino/StaMPS.git"
STAMPS_COMMIT = "7cabf05eddf8ebe8694e5346fe0f9d48aaef4962"
TRAIN_REPO = "https://github.com/pyccino/TRAIN.git"
TRAIN_COMMIT = "6d0273ae67d2a9f07a696b6a14298ef2c31607d8"
CORE_TOOLS = (
    "calamp", "cpxsum", "pscphase", "pscdem", "psclonlat",
    "selpsc_patch", "selsbc_patch",
)


def require_tools() -> None:
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "aarch64"}:
        raise RuntimeError("Linux StaMPS build requires x86_64 or aarch64 Linux.")
    missing = [name for name in ("git", "cmake", "ctest", "c++", "snaphu", "gawk", "csh")
               if not shutil.which(name)]
    if missing:
        raise RuntimeError("Install Linux prerequisites first: " + ", ".join(missing))


def fetch_pinned(repo: str, commit: str, destination: Path) -> None:
    subprocess.run(["git", "init", str(destination)], check=True)
    subprocess.run(["git", "-C", str(destination), "remote", "add", "origin", repo], check=True)
    subprocess.run(["git", "-C", str(destination), "fetch", "--depth", "1", "origin", commit],
                   check=True)
    subprocess.run(["git", "-C", str(destination), "checkout", "--detach", "FETCH_HEAD"],
                   check=True)
    actual = subprocess.check_output(["git", "-C", str(destination), "rev-parse", "HEAD"],
                                     text=True).strip()
    if actual != commit:
        raise RuntimeError(f"Unexpected third-party revision: {actual} (expected {commit}).")


def prepare(output: Path) -> Path:
    require_tools()
    output = output.expanduser().resolve()
    if output.exists():
        raise RuntimeError(f"Output already exists; refusing to replace it: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".phase-linux-runtime-", dir=output.parent) as tmp:
        stage = Path(tmp) / "runtime"
        stage.mkdir()
        stamps = stage / "StaMPS"
        train = stage / "TRAIN"
        fetch_pinned(STAMPS_REPO, STAMPS_COMMIT, stamps)
        fetch_pinned(TRAIN_REPO, TRAIN_COMMIT, train)
        build = Path(tmp) / "stamps-build"
        subprocess.run(["cmake", "-S", str(stamps / "src"), "-B", str(build),
                        "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_DISMPH=OFF"], check=True)
        subprocess.run(["cmake", "--build", str(build), "--parallel", "2"], check=True)
        subprocess.run(["ctest", "--test-dir", str(build), "--output-on-failure"], check=True)
        triangle_build = Path(tmp) / "triangle-build"
        subprocess.run(["cmake", "-S", str(stamps / "external" / "triangle"),
                        "-B", str(triangle_build), "-DCMAKE_BUILD_TYPE=Release"], check=True)
        subprocess.run(["cmake", "--build", str(triangle_build), "--parallel", "2"], check=True)
        for name in CORE_TOOLS:
            binary = stamps / "bin" / name
            if not binary.is_file() or not os.access(binary, os.X_OK):
                raise RuntimeError(f"StaMPS build did not produce {name}.")
        triangle = stamps / "external" / "triangle" / "bin" / "triangle"
        if not triangle.is_file() or not os.access(triangle, os.X_OK):
            raise RuntimeError("Triangle build did not produce an executable.")
        if not (train / "matlab" / "aps_linear.m").is_file():
            raise RuntimeError("TRAIN is missing aps_linear.m.")
        # Keep source and license notices, but not Git history in the installed runtime.
        shutil.rmtree(stamps / ".git")
        shutil.rmtree(train / ".git")
        (stage / "phase-runtime.json").write_text(json.dumps({
            "stampsCommit": STAMPS_COMMIT, "trainCommit": TRAIN_COMMIT,
            "platform": "linux-" + platform.machine()
        }, indent=2) + "\n", encoding="utf-8")
        stage.rename(output)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(prepare(args.output))


if __name__ == "__main__":
    main()
