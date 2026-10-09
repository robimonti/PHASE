#!/usr/bin/env python3
"""Prepare a local Apple Silicon StaMPS/TRAIN runtime for PHASE testing.

Use the pyccino/StaMPS fork used by PHASE on Windows. SNAPHU and Shewchuk's
Triangle and GNU awk must be obtained separately under their own terms; this script does
not redistribute or download them. Pass its StaMPS and TRAIN output folders
to install-phase-unix.py after the native build passes.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


CORE_TOOLS = (
    "calamp", "cpxsum", "pscphase", "pscdem", "psclonlat",
    "selpsc_patch", "selsbc_patch",
)


def arm64_executable(path: Path) -> bool:
    if not path.is_file() or not os.access(path, os.X_OK):
        return False
    result = subprocess.run(["lipo", "-archs", str(path)],
                            capture_output=True, text=True, check=False)
    return result.returncode == 0 and "arm64" in result.stdout.split()


def prepare(stamps_source: Path, train_source: Path, snaphu: Path,
            triangle: Path, gawk: Path, output: Path) -> Path:
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise RuntimeError("This runtime build requires an Apple Silicon Mac.")
    stamps_source = stamps_source.expanduser().resolve()
    train_source = train_source.expanduser().resolve()
    output = output.expanduser().resolve()
    if output.exists():
        raise RuntimeError(f"Output already exists; refusing to replace it: {output}")
    if any(output == root or root in output.parents or output in root.parents
           for root in (stamps_source, train_source)):
        raise RuntimeError("Output and source checkouts must be separate.")
    for required in ("src/CMakeLists.txt", "matlab/stamps.m", "bin/mt_prep_snap"):
        if not (stamps_source / required).is_file():
            raise RuntimeError(f"Incomplete StaMPS source: {required}")
    for required in ("matlab/aps_linear.m", "matlab/aps_weather_model.m",
                     "matlab/setparm_aps.m"):
        if not (train_source / required).is_file():
            raise RuntimeError(f"Incomplete TRAIN source: {required}")
    for label, binary in (("SNAPHU", snaphu), ("Triangle", triangle),
                          ("GNU awk", gawk)):
        if not arm64_executable(binary.expanduser().resolve()):
            raise RuntimeError(f"{label} must be an executable arm64 binary: {binary}")
    if not shutil.which("cmake") or not shutil.which("ctest"):
        raise RuntimeError("CMake and CTest are required to build StaMPS.")
    if not shutil.which("csh"):
        raise RuntimeError("StaMPS requires csh on the build/test Mac.")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".phase-macos-runtime-",
                                     dir=output.parent) as temporary:
        stage = Path(temporary) / "runtime"
        stage.mkdir()
        ignore = shutil.ignore_patterns(".git", "build", "*.o", ".DS_Store")
        stamps = stage / "StaMPS"
        train = stage / "TRAIN"
        shutil.copytree(stamps_source, stamps, ignore=ignore)
        shutil.copytree(train_source, train, ignore=ignore)
        build = Path(temporary) / "build"
        subprocess.run(["cmake", "-S", str(stamps / "src"), "-B", str(build),
                        "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_DISMPH=OFF"], check=True)
        subprocess.run(["cmake", "--build", str(build), "--config", "Release"],
                       check=True)
        test_env = os.environ.copy()
        test_env["PATH"] = str(gawk.expanduser().resolve().parent) + os.pathsep + test_env["PATH"]
        subprocess.run(["ctest", "--test-dir", str(build), "--output-on-failure"],
                       env=test_env, check=True)
        for name in CORE_TOOLS:
            binary = stamps / "bin" / name
            if not arm64_executable(binary):
                raise RuntimeError(f"StaMPS did not produce an arm64 {name} binary.")
        for source, destination in (
            (snaphu, stamps / "external" / "snaphu" / "bin" / "snaphu"),
            (triangle, stamps / "external" / "triangle" / "bin" / "triangle"),
            (gawk, stamps / "external" / "gawk" / "bin" / "gawk"),
        ):
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source.expanduser().resolve(), destination)
            destination.chmod(destination.stat().st_mode | 0o111)
        stage.rename(output)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stamps-source", required=True, type=Path)
    parser.add_argument("--train-source", required=True, type=Path)
    parser.add_argument("--snaphu", required=True, type=Path)
    parser.add_argument("--triangle", required=True, type=Path)
    parser.add_argument("--gawk", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(prepare(args.stamps_source, args.train_source, args.snaphu,
                  args.triangle, args.gawk, args.output))


if __name__ == "__main__":
    main()
