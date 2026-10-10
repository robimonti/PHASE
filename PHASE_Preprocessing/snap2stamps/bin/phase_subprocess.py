"""Silent, live-output subprocess helpers for PHASE SNAP scripts."""

from __future__ import print_function

import os
from pathlib import Path
import re
import subprocess
import sys


_ORIGINAL_POPEN = subprocess.Popen


def live_popen(*args, **kwargs):
    """Create a Popen-compatible process whose communicate() tees output.

    On Windows CREATE_NO_WINDOW prevents SNAP gpt.exe from flashing separate
    console windows.  Existing PHASE scripts can keep using communicate()[0]
    and returncode while their output is streamed to MATLAB at the same time.
    """

    if sys.platform == "darwin" and args and isinstance(args[0], (list, tuple)):
        command, notice = macos_gpt_command(args[0])
        args = (command,) + args[1:]
        if notice:
            print(notice, flush=True)
    if os.name == "nt":
        kwargs["creationflags"] = (
            int(kwargs.get("creationflags", 0))
            | int(getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000))
        )
    return _LiveProcess(_ORIGINAL_POPEN(*args, **kwargs))


def _size_bytes(value):
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([KMG])\s*", str(value), re.I)
    if not match:
        return None
    return int(float(match.group(1)) * 1024 ** {"K": 1, "M": 2, "G": 3}[match.group(2).upper()])


def _mac_physical_memory():
    try:
        return int(subprocess.check_output(
            ["/usr/sbin/sysctl", "-n", "hw.memsize"], text=True, timeout=3).strip())
    except (OSError, ValueError, subprocess.SubprocessError):
        return 8 * 1024 ** 3  # Conservative when macOS cannot report RAM.


def _gpt_heap_bytes(gpt):
    options = Path(str(gpt)).with_name("gpt.vmoptions")
    try:
        content = options.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return 5 * 1024 ** 3
    values = re.findall(r"(?m)^\s*-Xmx\s*(\d+(?:\.\d+)?\s*[KMG])\s*$", content, re.I)
    return (_size_bytes(values[-1]) if values else None) or 5 * 1024 ** 3


def macos_gpt_command(command, physical_bytes=None, heap_bytes=None):
    """Keep tile cache and concurrency below a Mac's real SNAP memory budget.

    The PHASE project configuration is not rewritten; the effective values are
    logged for each GPT invocation. This changes performance, not processing
    operators or scientific parameters.
    """
    adjusted = list(command)
    if not adjusted or Path(str(adjusted[0])).name.lower() != "gpt":
        return adjusted, ""
    physical = physical_bytes or _mac_physical_memory()
    heap = heap_bytes or _gpt_heap_bytes(adjusted[0])
    mebibyte = 1024 ** 2
    budget = min(int(physical * 0.08), int(heap * 0.20))
    safe_mebibytes = max(128, 2 ** max(7, (budget // mebibyte).bit_length() - 1))
    safe_cache = safe_mebibytes * mebibyte
    safe_threads = 2 if physical <= 12 * 1024 ** 3 else 4 if physical <= 24 * 1024 ** 3 else 8
    changes = []
    for flag in ("-c", "-q"):
        if flag not in adjusted:
            continue
        index = adjusted.index(flag) + 1
        if index >= len(adjusted):
            continue
        original = str(adjusted[index])
        if flag == "-c":
            requested = _size_bytes(original)
            if requested is None or requested > safe_cache:
                adjusted[index] = "%dM" % safe_mebibytes
        else:
            try:
                requested = int(original)
            except ValueError:
                requested = safe_threads + 1
            if requested > safe_threads:
                adjusted[index] = str(safe_threads)
        if str(adjusted[index]) != original:
            changes.append("%s %s → %s" % (flag, original, adjusted[index]))
    if "-x" not in adjusted:
        adjusted.append("-x")
        changes.append("-x (clear completed tile rows)")
    notice = "PHASE macOS SNAP memory guard: " + ", ".join(changes) if changes else ""
    return adjusted, notice


class _LiveProcess(object):
    def __init__(self, process):
        self._process = process

    @property
    def returncode(self):
        return self._process.returncode

    @property
    def stdout(self):
        return self._process.stdout

    def communicate(self, input=None, timeout=None):
        if input is not None or timeout is not None or self._process.stdout is None:
            return self._process.communicate(input=input, timeout=timeout)

        chunks = []
        while True:
            chunk = self._process.stdout.readline()
            if not chunk:
                if self._process.poll() is not None:
                    break
                continue
            chunks.append(chunk)
            _forward(chunk)
        self._process.wait()
        return b"".join(chunks), None

    def __getattr__(self, name):
        return getattr(self._process, name)


def _forward(chunk):
    stream = getattr(sys.stdout, "buffer", None)
    if stream is not None:
        stream.write(chunk)
        stream.flush()
        return
    if sys.stdout is not None:
        sys.stdout.write(chunk.decode("utf-8", errors="replace"))
        sys.stdout.flush()
