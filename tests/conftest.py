"""Shared pytest fixtures and marker-driven skip logic for the PHASE suite."""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest


@pytest.fixture
def phase_root() -> Path:
    """Absolute path to the PHASE repository root."""
    return Path(__file__).parent.parent


def pytest_collection_modifyitems(config: pytest.Config,
                                  items: list[pytest.Item]) -> None:
    """Auto-skip tests marked windows_only on non-Windows and requires_matlab
    when the `matlab` executable is not on PATH."""
    skip_non_windows = pytest.mark.skip(
        reason="test marked windows_only; current platform is not Windows")
    skip_no_matlab = pytest.mark.skip(
        reason="test marked requires_matlab; `matlab` not found on PATH")
    have_matlab = shutil.which("matlab") is not None
    is_windows = sys.platform.startswith("win")
    for item in items:
        if "windows_only" in item.keywords and not is_windows:
            item.add_marker(skip_non_windows)
        if "requires_matlab" in item.keywords and not have_matlab:
            item.add_marker(skip_no_matlab)
