"""Silent, live-output subprocess helpers for PHASE SNAP scripts."""

from __future__ import print_function

import os
import subprocess
import sys


_ORIGINAL_POPEN = subprocess.Popen


def live_popen(*args, **kwargs):
    """Create a Popen-compatible process whose communicate() tees output.

    On Windows CREATE_NO_WINDOW prevents SNAP gpt.exe from flashing separate
    console windows.  Existing PHASE scripts can keep using communicate()[0]
    and returncode while their output is streamed to MATLAB at the same time.
    """

    if os.name == "nt":
        kwargs["creationflags"] = (
            int(kwargs.get("creationflags", 0))
            | int(getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000))
        )
    return _LiveProcess(_ORIGINAL_POPEN(*args, **kwargs))


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
