"""Detect SNAP-version mismatch between BEAM-DIMAP products and the gpt
binary about to read them.

Why this exists
---------------
Starting with SNAP 13, the StampsExport operator chokes on `.dim` products
written by SNAP 9 with a NullPointerException in
`StampsExportOp.initialize`:

    Cannot invoke "ProductData.getElems()" because the return value of
    "TiePointGrid.getGridData()" is null

Root cause: the tie-point grid serialization changed between SNAP 9 and
SNAP 13, and SNAP 13's `InputProductValidator.checkIfCompatibleProducts`
silently fails to load the legacy TPG, then trips on the null. The
operator returns no clear error to the user; downstream StaMPS hangs
mid-PSI without a usable diagnostic.

There is no upstream fix. The only reliable workaround is to re-run the
preprocessing from scratch with whichever SNAP version the user has
installed. This module flags the mismatch up-front.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

# What SNAP actually writes in BEAM-DIMAP `.dim` files
# ------------------------------------------------------
# Each Processing_Graph node carries a `moduleVersion` MDATTR which tracks
# the version of the SNAP framework / S1TBX module that ran the step:
#   <MDATTR name="moduleVersion" type="ascii" mode="rw">9.0.0</MDATTR>
#   <MDATTR name="moduleVersion" type="ascii" mode="rw">13.0.0</MDATTR>
# All nodes of a single `.dim` carry the same value (= the SNAP build that
# produced the product). We pick the highest moduleVersion seen to be
# robust against products that mixed third-party operators with older
# embedded versions.
#
# Two earlier regex hypotheses (`processing_software_version`,
# `snap_version`) do NOT match real SNAP output — those keys exist in some
# documentation examples but not in actual `.dim` files emitted by gpt.
_DIM_MODULE_VERSION_RX = re.compile(
    r'<MDATTR\s+name="moduleVersion"[^>]*>\s*'
    r'(?P<ver>\d+\.\d+(?:\.\d+)?)\s*</MDATTR>',
    re.IGNORECASE,
)

# Legacy fallback: kept so synthetic test fixtures and any hand-edited
# `.dim` using these keys still get recognized.
_DIM_VERSION_LEGACY_RX = re.compile(
    r'<MDATTR\s+name="processing_software_version"[^>]*>\s*'
    r'(?P<ver>\d+\.\d+(?:\.\d+)?)\s*</MDATTR>',
    re.IGNORECASE,
)
_DIM_VERSION_LEGACY_RX_ALT = re.compile(
    r'<MDATTR[^>]*name="snap_version"[^>]*>\s*(?P<ver>\d+\.\d+(?:\.\d+)?)',
    re.IGNORECASE,
)


def _version_tuple(ver: str) -> tuple[int, ...]:
    return tuple(int(x) for x in ver.split("."))


def detect_dim_snap_version(dim_path: Path) -> str | None:
    """Return the SNAP version that wrote the `.dim` file, or None if
    not detectable. Reads only the XML header — no `.data/` access."""
    if not dim_path.is_file():
        return None
    try:
        text = dim_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    versions = [m.group("ver") for m in _DIM_MODULE_VERSION_RX.finditer(text)]
    if versions:
        return max(versions, key=_version_tuple)
    for rx in (_DIM_VERSION_LEGACY_RX, _DIM_VERSION_LEGACY_RX_ALT):
        m = rx.search(text)
        if m:
            return m.group("ver")
    return None


def detect_gpt_snap_version(gpt_path: Path) -> str | None:
    """Run `gpt --diag` and parse the `SNAP Release version X.Y.Z` line."""
    if not gpt_path.is_file():
        return None
    try:
        result = subprocess.run(
            [str(gpt_path), "--diag"],
            capture_output=True, text=True, timeout=60, check=False,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None
    blob = (result.stdout or "") + (result.stderr or "")
    m = re.search(
        r"SNAP Release version\s+(?P<ver>\d+\.\d+(?:\.\d+)?)", blob)
    return m.group("ver") if m else None


def major(version: str) -> int:
    """First dotted component as int. '13.0.0' -> 13."""
    return int(version.split(".")[0])


def check_compatibility(
    dim_paths: list[Path],
    gpt_path: Path,
) -> dict[str, list[str]]:
    """Cross-check each `.dim` against the gpt SNAP version.

    Returns:
        {'matched': [...], 'mismatched': [(dim, dim_ver, gpt_ver)],
         'unknown': [dim]}
    """
    gpt_ver = detect_gpt_snap_version(gpt_path)
    report: dict[str, list] = {
        "matched": [], "mismatched": [], "unknown": [], "gpt_version": gpt_ver,
    }
    if gpt_ver is None:
        # Cannot determine gpt version; skip the check rather than block.
        return report
    for dim in dim_paths:
        dim_ver = detect_dim_snap_version(dim)
        if dim_ver is None:
            report["unknown"].append(str(dim))
            continue
        if major(dim_ver) == major(gpt_ver):
            report["matched"].append(str(dim))
        else:
            report["mismatched"].append((str(dim), dim_ver, gpt_ver))
    return report


def warn_if_mismatch(report: dict, *, write=print) -> bool:
    """Print a human-actionable warning if mismatches are present.
    Returns True iff any mismatch was found."""
    mismatches = report.get("mismatched", [])
    if not mismatches:
        return False
    gpt_ver = report.get("gpt_version", "?")
    write("=" * 70)
    write(f"WARN: SNAP version mismatch detected")
    write(f"      gpt is SNAP {gpt_ver}, but {len(mismatches)} .dim "
          f"file(s) were produced by an older SNAP:")
    for dim, dim_ver, _ in mismatches:
        write(f"        {dim_ver}: {dim}")
    write("")
    write("      SNAP 13's StampsExport NullPointerExceptions on legacy")
    write("      tie-point grids from SNAP 9 .dim products. StaMPS will")
    write("      then hang mid-PSI without a clear error.")
    write("")
    write("      Recommended fix: re-run the preprocessing pipeline")
    write("      (split + coreg + ifg) from the original .SAFE.zip with")
    write("      the same SNAP version that will run StampsExport. Do not")
    write("      mix .dim files from different SNAP majors.")
    write("=" * 70)
    return True
