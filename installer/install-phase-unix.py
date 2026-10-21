#!/usr/bin/env python3
"""Per-user PHASE Hub installer for Apple Silicon macOS and Linux.

MATLAB and SNAP are proprietary/external prerequisites. StaMPS and TRAIN may
be supplied as already prepared runtimes; their scientific execution is not
claimed to be validated by this installer.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import plistlib
import re
import shutil
import shlex
import subprocess
import sys
import tempfile
import venv


PHASE_REPO = "https://github.com/robimonti/PHASE.git"
STAMPS_COMMIT = "7cabf05eddf8ebe8694e5346fe0f9d48aaef4962"
TRAIN_COMMIT = "6d0273ae67d2a9f07a696b6a14298ef2c31607d8"
PYTHON_PACKAGES = ("openpyxl", "requests", "asf_search", "shapely", "certifi")
MACOS_PSI_TOOLS = (
    "bin/calamp", "bin/cpxsum", "bin/pscphase", "bin/pscdem",
    "bin/psclonlat", "bin/selpsc_patch", "bin/selsbc_patch",
    "external/snaphu/bin/snaphu", "external/triangle/bin/triangle",
    "external/gawk/bin/gawk",
)
REQUIRED_FILES = (
    "PHASE_Hub.m",
    "PHASE_Hub_UI.html",
    "PHASE_Stamps_Dataset_UI.html",
    "PHASE_logo.png",
    "PHASE_mod1a.png",
    "PHASE_mod1b.png",
    "PHASE_mod2.png",
    "+phase_hub/App.m",
    "+phase_project/open.m",
    "+phase_project/displayVersion.m",
    "+phase_project/workflowStatus.m",
    "phase_update.py",
    "PHASE_Preprocessing_beta.m",
    "PHASE_Preprocessing/+phase_preprocessing_beta/App.m",
    "PHASE_Preprocessing/+phase_preprocessing_beta/ProcessingEngine.m",
    "PHASE_Preprocessing/+phase_stamps_beta/App.m",
    "PHASE_Preprocessing/+phase_stamps_beta/ProcessingContext.m",
    "+phase_model_beta/App.m",
    "+phase_model_beta/ProcessingEngine.m",
)
ROOT_IGNORED = {
    ".git", ".github", ".gitignore", ".pytest_cache", ".codex", ".agents", "legacy",
    "tests", "installer", "pytest.ini", "prepare-windows-runtime.ps1"
}


def system_name() -> str:
    if sys.platform == "darwin":
        if platform.machine() != "arm64":
            raise RuntimeError("PHASE for macOS requires Apple Silicon (arm64).")
        return "macos"
    if sys.platform.startswith("linux"):
        return "linux"
    raise RuntimeError("Use install-phase.ps1 on Windows.")


def default_prefix(system: str) -> Path:
    if system == "macos":
        return Path.home() / "Library" / "Application Support" / "PHASE"
    return Path.home() / ".local" / "share" / "PHASE"


def find_matlab(system: str, specified: str | None) -> Path:
    if specified:
        candidates = [Path(specified).expanduser()]
    elif system == "macos":
        candidates = sorted(Path("/Applications").glob("MATLAB_R*.app/bin/matlab"), reverse=True)
    else:
        path_binary = shutil.which("matlab")
        candidates = [Path(path_binary)] if path_binary else []
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise RuntimeError("MATLAB not found. Supply --matlab /path/to/matlab.")


def find_gpt(system: str, specified: str | None) -> Path:
    if specified:
        candidates = [Path(specified).expanduser()]
    elif system == "macos":
        candidates = [Path("/Applications/esa-snap/bin/gpt")]
    else:
        candidates = [Path("/opt/esa-snap/bin/gpt"), Path("/opt/snap/bin/gpt")]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK) and candidate.name == "gpt":
            return candidate.resolve()
    raise RuntimeError(
        "ESA SNAP gpt not found. Install SNAP and supply --gpt /path/to/snap/bin/gpt "
        "(the unrelated /usr/sbin/gpt is not SNAP)."
    )


def find_python(specified: str | None) -> Path:
    candidate = Path(specified).expanduser() if specified else Path(sys.executable)
    if not candidate.is_file():
        raise RuntimeError(f"Python executable not found: {candidate}")
    output = subprocess.check_output(
        [str(candidate), "-c", "import sys; print(sys.version_info[0], sys.version_info[1])"],
        text=True,
    ).strip()
    major, minor = map(int, output.split())
    if major != 3 or minor < 10:
        raise RuntimeError("PHASE requires Python 3.10 or newer.")
    return candidate.resolve()


def validate_runtime(source: Path) -> None:
    missing = [name for name in REQUIRED_FILES if not (source / name).is_file()]
    if missing:
        raise RuntimeError("Incomplete PHASE source: " + ", ".join(missing))


def validate_external(path: str | None, kind: str) -> Path | None:
    if not path:
        return None
    root = Path(path).expanduser().resolve()
    required = ("matlab/stamps.m", "matlab/setparm.m") if kind == "stamps" else ("matlab",)
    if any(not (root / name).exists() for name in required):
        raise RuntimeError(f"Incomplete {kind} runtime: {root}")
    if kind == "stamps":
        triangle = root / "external" / "triangle" / "bin" / "triangle"
        if not triangle.is_file() or not os.access(triangle, os.X_OK):
            raise RuntimeError(f"StaMPS Triangle is not executable: {triangle}")
        result = subprocess.run([str(triangle), "-h"], capture_output=True,
                                text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"StaMPS Triangle cannot start: {triangle}")
    return root


def macos_psi_missing(stamps: Path | None, train: Path | None) -> list[str]:
    """Report prerequisites before claiming an Apple Silicon PSI runtime."""
    def native(path: Path) -> bool:
        if not path.is_file() or not os.access(path, os.X_OK):
            return False
        result = subprocess.run(["lipo", "-archs", str(path)],
                                capture_output=True, text=True, check=False)
        return result.returncode == 0 and "arm64" in result.stdout.split()

    missing = []
    if stamps is None:
        missing.append("StaMPS runtime")
    else:
        for name in MACOS_PSI_TOOLS:
            if not native(stamps / name):
                missing.append(f"StaMPS/{name}")
    if train is None or not (train / "matlab" / "aps_linear.m").is_file():
        missing.append("TRAIN runtime (optional if correction is disabled)")
    return missing


def trust_bundled_macos_tools(stamps_source: Path, train_source: Path,
                              staged_stamps: Path) -> int:
    """Clear quarantine only for explicitly approved, pinned PSI executables.

    This is an opt-in alternative to macOS's per-tool Open Anyway prompts for
    the unsigned release. It does not modify Gatekeeper or unrelated files.
    """
    if stamps_source.parent != train_source.parent:
        raise RuntimeError("StaMPS and TRAIN must belong to the same prepared runtime.")
    manifest = stamps_source.parent / "phase-runtime.json"
    if not manifest.is_file():
        raise RuntimeError("Trusted-tool option requires a pinned PHASE runtime manifest.")
    info = json.loads(manifest.read_text(encoding="utf-8"))
    if (info.get("platform") != "macos-arm64" or
            info.get("stampsCommit") != STAMPS_COMMIT or
            info.get("trainCommit") != TRAIN_COMMIT):
        raise RuntimeError("Trusted-tool option requires the pinned macOS PSI runtime.")
    hashes = info.get("binarySha256")
    if not isinstance(hashes, dict) or set(hashes) != set(MACOS_PSI_TOOLS):
        raise RuntimeError("Trusted-tool option requires checksums for every PSI executable.")
    cleared = 0
    for relative in MACOS_PSI_TOOLS:
        binary = staged_stamps / relative
        if binary.is_symlink() or not binary.is_file() or not os.access(binary, os.X_OK):
            raise RuntimeError(f"Invalid bundled PSI executable: {relative}")
        digest = hashlib.sha256(binary.read_bytes()).hexdigest()
        if digest != hashes[relative]:
            raise RuntimeError(f"Bundled PSI executable checksum mismatch: {relative}")
        arches = subprocess.run(["lipo", "-archs", str(binary)],
                                capture_output=True, text=True, check=True)
        if "arm64" not in arches.stdout.split():
            raise RuntimeError(f"Not an Apple Silicon PSI executable: {relative}")
        attributes = subprocess.run(["/usr/bin/xattr", str(binary)],
                                    capture_output=True, text=True, check=True)
        if "com.apple.quarantine" in attributes.stdout.splitlines():
            subprocess.run(["/usr/bin/xattr", "-d", "com.apple.quarantine",
                            str(binary)], check=True)
            cleared += 1
    return cleared


def ignored(directory: str, names: list[str]) -> set[str]:
    current = Path(directory)
    blocked = {name for name in names if name in {"__pycache__", ".DS_Store"}
               or name.endswith((".pyc", ".docx", ".zip"))}
    if current == SOURCE_FOR_COPY:
        blocked.update(ROOT_IGNORED.intersection(names))
    if current.name == "tools":
        blocked.update(name for name in names if name not in {
            "srtm_precache.py", "snap_dim_version_check.py"
        })
    return blocked


SOURCE_FOR_COPY = Path("/")


def write_launcher(path: Path, prefix: Path, matlab: Path, gpt: Path, python: Path) -> None:
    engine = prefix / "engine"
    stamps = engine / "StaMPS"
    train = engine / "TRAIN"
    matlab_root = str(engine).replace("'", "''")
    stamps_matlab = str(stamps / "matlab").replace("'", "''")
    train_root = str(train).replace("'", "''")
    train_matlab = str(train / "matlab").replace("'", "''")
    expression = (
        f"addpath('{matlab_root}'); "
        f"addpath(fullfile('{matlab_root}','PHASE_Preprocessing')); "
        f"if isfolder('{stamps_matlab}'), addpath('{stamps_matlab}'); end; "
        f"if isfolder('{train_matlab}'), addpath('{train_root}'); "
        f"addpath(genpath('{train_matlab}')); end; "
        "PHASE_Hub"
    )
    content = (
        "#!/bin/sh\nset -eu\n"
        f"export PHASE_GPTBIN={shlex.quote(str(gpt))}\n"
        f"export PHASE_PYTHON={shlex.quote(str(python))}\n"
        f"export STAMPS={shlex.quote(str(stamps))}\n"
        f"export APS_toolbox={shlex.quote(str(train))}\n"
        f"export PATH=/opt/homebrew/bin:{shlex.quote(str(matlab.parent))}:"
        f"{shlex.quote(str(stamps / 'bin'))}:"
        f"{shlex.quote(str(stamps / 'external' / 'snaphu' / 'bin'))}:"
        f"{shlex.quote(str(stamps / 'external' / 'triangle' / 'bin'))}:"
        f"{shlex.quote(str(stamps / 'external' / 'gawk' / 'bin'))}:"
        f"{shlex.quote(str(python.parent))}:\"$PATH\"\n"
        f"{shlex.quote(str(python))} {shlex.quote(str(engine / 'phase_update.py'))} "
        f"apply --prefix {shlex.quote(str(prefix))}\n"
        f"exec {shlex.quote(str(matlab))} -desktop -r {shlex.quote(expression)}\n"
    )
    path.write_text(content, encoding="utf-8")
    path.chmod(0o755)


def configure_unix_runtimes(stage: Path, prefix: Path) -> None:
    """Replace upstream example paths with relocatable, per-user paths."""
    stamps = stage / "engine" / "StaMPS"
    if stamps.is_dir():
        (stamps / "StaMPS_CONFIG.bash").write_text(
            "# PHASE-managed StaMPS environment.\n"
            f"export STAMPS={shlex.quote(str(prefix / 'engine' / 'StaMPS'))}\n"
            'export MATLABPATH="$STAMPS/matlab${MATLABPATH:+:$MATLABPATH}"\n'
            'export PYTHONPATH="$STAMPS/python${PYTHONPATH:+:$PYTHONPATH}"\n'
            'export PATH="$STAMPS/bin:$STAMPS/external/snaphu/bin:'
            '$STAMPS/external/triangle/bin:$STAMPS/external/gawk/bin:$PATH"\n',
            encoding="utf-8",
        )
    train = stage / "engine" / "TRAIN"
    if train.is_dir():
        (train / "APS_CONFIG.sh").write_text(
            "# PHASE-managed TRAIN environment.\n"
            f"export APS_toolbox={shlex.quote(str(prefix / 'engine' / 'TRAIN'))}\n"
            'export MATLABPATH="$APS_toolbox/matlab${MATLABPATH:+:$MATLABPATH}"\n'
            'export PATH="$APS_toolbox/bin:$PATH"\n',
            encoding="utf-8",
        )


def release_tag(checkout: Path) -> str:
    """Return an exact v7 release tag, or an empty string for development builds."""
    result = subprocess.run(["git", "describe", "--tags", "--exact-match"],
                            cwd=checkout, capture_output=True, text=True, check=False)
    tag = result.stdout.strip() if result.returncode == 0 else ""
    return tag if re.fullmatch(r"v7\.\d+\.\d+", tag) else ""


def embed_release_metadata(engine: Path, checkout: Path) -> None:
    """Mark an installer payload only when built from an exact v7 Git tag."""
    tag = release_tag(checkout)
    if tag:
        (engine / "phase-release.json").write_text(json.dumps({
            "tag": tag, "updateSchema": 1
        }, indent=2) + "\n", encoding="utf-8")


def payload_version(engine: Path) -> str:
    manifest = engine / "phase-release.json"
    if not manifest.is_file():
        return "dev"
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
        tag = data.get("tag", "")
        if data.get("updateSchema") == 1 and re.fullmatch(r"v7\.\d+\.\d+", tag):
            return tag
    except (OSError, ValueError, TypeError):
        pass
    raise RuntimeError("Invalid PHASE release metadata in installer payload.")


def create_macos_app(stage: Path, prefix: Path) -> None:
    bundle = stage / "PHASE.app" / "Contents"
    executable = bundle / "MacOS" / "PHASE"
    executable.parent.mkdir(parents=True)
    executable.write_text(
        "#!/bin/sh\n" + f"exec {shlex.quote(str(prefix / 'launch-phase.sh'))} \"$@\"\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    icon = Path(__file__).resolve().parent / "PHASE.icns"
    if icon.is_file():
        resources = bundle / "Resources"
        resources.mkdir()
        shutil.copy2(icon, resources / "PHASE.icns")
    with (bundle / "Info.plist").open("wb") as stream:
        version = payload_version(stage / "engine")
        if version == "dev":
            version = "v7.0.2"
        plistlib.dump({
            "CFBundleName": "PHASE", "CFBundleDisplayName": "PHASE",
            "CFBundleIdentifier": "org.phaseinsar.phase", "CFBundleVersion": version[1:],
            "CFBundleShortVersionString": version[1:], "CFBundleExecutable": "PHASE",
            "CFBundlePackageType": "APPL", "LSMinimumSystemVersion": "13.0",
            "CFBundleIconFile": "PHASE.icns",
        }, stream)


def create_linux_desktop(stage: Path, prefix: Path) -> None:
    executable = str(prefix / "launch-phase.sh").replace('\\', '\\\\').replace('"', '\\"')
    (stage / "PHASE.desktop").write_text(
        "[Desktop Entry]\nType=Application\nName=PHASE\n"
        f'Exec="{executable}"\nIcon={prefix / "engine" / "Logo_square.png"}\n'
        "Terminal=false\nCategories=Science;Education;\n",
        encoding="utf-8",
    )


def ensure_safe_prefix(prefix: Path, source: Path | None) -> None:
    prohibited = {Path("/"), Path.home().resolve()}
    if prefix in prohibited:
        raise RuntimeError(f"Unsafe installation prefix: {prefix}")
    if source and (prefix == source or prefix in source.parents or source in prefix.parents):
        raise RuntimeError("Installation prefix and source checkout must be separate.")
    if prefix.is_dir() and any(prefix.iterdir()) and not (prefix / "install.json").is_file():
        raise RuntimeError(
            f"{prefix} contains files but is not marked as a managed PHASE installation."
        )


def create_shortcut(system: str, prefix: Path) -> str:
    if system == "macos":
        folder = Path.home() / "Applications"
        link = folder / "PHASE.app"
        target = prefix / "PHASE.app"
    else:
        folder = Path.home() / ".local" / "share" / "applications"
        link = folder / "phase.desktop"
        target = prefix / "PHASE.desktop"
    folder.mkdir(parents=True, exist_ok=True)
    if link.is_symlink() and link.resolve() == target:
        return str(link)
    if link.exists() or link.is_symlink():
        return f"Shortcut not replaced (existing unrelated item): {link}; launch {target} manually."
    link.symlink_to(target)
    if system == "linux":
        bin_folder = Path.home() / ".local" / "bin"
        bin_folder.mkdir(parents=True, exist_ok=True)
        command = bin_folder / "phase"
        if not command.exists() and not command.is_symlink():
            command.symlink_to(prefix / "launch-phase.sh")
    return str(link)


def install(args: argparse.Namespace) -> dict[str, str]:
    system = system_name()
    prefix = Path(args.prefix).expanduser().resolve() if args.prefix else default_prefix(system)
    matlab = find_matlab(system, args.matlab)
    gpt = find_gpt(system, args.gpt)
    python = find_python(args.python)
    source = Path(args.source).expanduser().resolve() if args.source else None
    if source:
        validate_runtime(source)
    stamps = validate_external(args.stamps, "stamps")
    train = validate_external(args.train, "train")
    if args.trust_bundled_macos_tools and (system != "macos" or not stamps or not train):
        raise RuntimeError("Trusted-tool option requires the bundled macOS StaMPS/TRAIN runtime.")
    psi_missing = macos_psi_missing(stamps, train) if system == "macos" else []
    if stamps and not args.dry_run and any(item != "TRAIN runtime (optional if correction is disabled)"
                                      for item in psi_missing):
        raise RuntimeError("Incomplete macOS PSI runtime: " + ", ".join(psi_missing))
    ensure_safe_prefix(prefix, source)
    plan = {
        "system": system, "prefix": str(prefix),
        "source": str(source) if source else f"{PHASE_REPO}#{args.branch}",
        "matlab": str(matlab), "gpt": str(gpt), "python": str(python),
        "stamps": str(stamps) if stamps else "not supplied",
        "train": str(train) if train else "not supplied",
        "psiRuntimeMissing": psi_missing,
        "trustedBundledTools": bool(args.trust_bundled_macos_tools),
        "pythonDependencies": "skipped" if args.skip_python_deps else ", ".join(PYTHON_PACKAGES),
    }
    if args.dry_run:
        return plan

    prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".phase-install-", dir=prefix.parent) as temporary:
        stage = Path(temporary) / "payload"
        stage.mkdir()
        if source is None:
            source = Path(temporary) / "source"
            subprocess.run([
                "git", "clone", "--depth", "1", "--branch", args.branch,
                PHASE_REPO, str(source)
            ], check=True)
            validate_runtime(source)
        global SOURCE_FOR_COPY
        SOURCE_FOR_COPY = source
        shutil.copytree(source, stage / "engine", ignore=ignored)
        if stamps:
            shutil.copytree(stamps, stage / "engine" / "StaMPS")
        if train:
            shutil.copytree(train, stage / "engine" / "TRAIN")
        if args.trust_bundled_macos_tools:
            plan["quarantinedToolsApproved"] = str(trust_bundled_macos_tools(
                stamps, train, stage / "engine" / "StaMPS"))
        configure_unix_runtimes(stage, prefix)
        if args.skip_python_deps:
            installed_python = python
        else:
            venv.EnvBuilder(with_pip=True).create(stage / "venv")
            staged_python = stage / "venv" / "bin" / "python"
            subprocess.run([
                str(staged_python), "-m", "pip", "install", "--disable-pip-version-check",
                "--no-input", *PYTHON_PACKAGES
            ], check=True)
            installed_python = prefix / "venv" / "bin" / "python"
        write_launcher(stage / "launch-phase.sh", prefix, matlab, gpt, installed_python)
        if system == "macos":
            create_macos_app(stage, prefix)
        else:
            create_linux_desktop(stage, prefix)
        (stage / "install.json").write_text(json.dumps({
            **plan, "installedAt": datetime.now(timezone.utc).isoformat(),
            "updateSchema": 1, "version": payload_version(stage / "engine")
        }, indent=2) + "\n", encoding="utf-8")

        prefix.mkdir(parents=True, exist_ok=True)
        backup = prefix / "backups" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        previous: list[tuple[Path, Path]] = []
        added: list[Path] = []
        try:
            for item in stage.iterdir():
                target = prefix / item.name
                if target.exists():
                    backup.mkdir(parents=True, exist_ok=True)
                    saved = backup / item.name
                    if saved.exists():
                        raise RuntimeError(f"Backup collision: {saved}")
                    target.rename(saved)
                    previous.append((target, saved))
                item.rename(target)
                added.append(target)
        except Exception:
            for target in reversed(added):
                if target.is_dir():
                    shutil.rmtree(target)
                elif target.exists():
                    target.unlink()
            for target, saved in reversed(previous):
                saved.rename(target)
            raise
    if args.no_shortcut:
        plan["shortcut"] = "skipped"
    else:
        try:
            plan["shortcut"] = create_shortcut(system, prefix)
        except OSError as error:
            plan["shortcut"] = f"Could not create shortcut: {error}"
    plan["launcher"] = str(prefix / "launch-phase.sh")
    return plan


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", help="Local PHASE checkout; default clones GitHub")
    parser.add_argument("--branch", default="main", help="Git branch for remote installation")
    parser.add_argument("--prefix", help="Per-user installation directory")
    parser.add_argument("--matlab", help="MATLAB executable")
    parser.add_argument("--gpt", help="ESA SNAP gpt executable")
    parser.add_argument("--python", help="Python 3.10+ executable")
    parser.add_argument("--stamps", help="Prepared StaMPS runtime to bundle")
    parser.add_argument("--train", help="Prepared TRAIN runtime to bundle")
    parser.add_argument("--trust-bundled-macos-tools", action="store_true",
                        help="With explicit user consent, clear quarantine only on the ten pinned bundled PSI executables")
    parser.add_argument("--skip-python-deps", action="store_true", help="Development only")
    parser.add_argument("--no-shortcut", action="store_true", help="Do not create an app-menu shortcut")
    parser.add_argument("--dry-run", action="store_true", help="Validate without writing")
    args = parser.parse_args()
    try:
        print(json.dumps(install(args), indent=2, ensure_ascii=False))
        return 0
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"PHASE installation failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
