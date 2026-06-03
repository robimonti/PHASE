"""Pre-cache SRTM 3Sec DEM tiles in SNAP's auxdata folder before running
StampsExport.

Why this exists
---------------
SNAP's StampsExport operator has SRTM 3Sec hardcoded for generating the
geocoding lat/lon files, regardless of the DEM the user chose for
coregistration. When gpt runs offline (or the auto-download fails for any
reason — SNAP 13 is particularly fragile here on Sentinel-1C/1D inputs),
StampsExport silently produces partially corrupted geo files. StaMPS then
hangs mid-PSI without a clear error.

This pre-cache step downloads the SRTM 3Sec tiles covering the project AOI
into SNAP's auxdata folder *before* gpt runs, eliminating the auto-download
path that breaks.

CLI:
    python tools/srtm_precache.py <project.conf>

Reads LATMIN/LATMAX/LONMIN/LONMAX from the conf, computes the covering
SRTM 3Sec tiles, downloads any missing ones into:
    %USERPROFILE%/.snap/auxdata/dem/SRTM 3Sec/   (Windows)
    ~/.snap/auxdata/dem/SRTM 3Sec/               (Unix)
"""
from __future__ import annotations

import math
import os
import sys
import urllib.request
import urllib.error
from pathlib import Path

# CGIAR-CSI SRTM 3Sec is distributed as 5°x5° tiles indexed by column/row.
# Origin: top-left (lon=-180°, lat=+60°). Tile size: 5°.
# Naming: srtm_<COL>_<ROW>.zip where COL: 1..72, ROW: 1..24.
SRTM_BASE_URL = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF"

# Mirror fallback if primary fails
SRTM_MIRRORS = [
    SRTM_BASE_URL,
    "https://download.osgeo.org/dem/srtm/srtm_5x5/TIFF",  # OSGeo mirror (community)
]


def latlon_to_tile(lat: float, lon: float) -> tuple[int, int]:
    """Map a (lat, lon) point to its SRTM 3Sec tile (col, row).

    Origin is the top-left corner: lon=-180°, lat=+60°. Tiles are 5° wide.
    The tile that *contains* a point at the bottom or right edge belongs
    to the next index by SNAP/CGIAR convention; we use floor() which
    matches that.
    """
    col = math.floor((lon + 180.0) / 5.0) + 1
    row = math.floor((60.0 - lat) / 5.0) + 1
    return col, row


def aoi_tiles(latmin: float, latmax: float, lonmin: float,
              lonmax: float) -> list[tuple[int, int]]:
    """Return the set of (col, row) SRTM 3Sec tiles covering the AOI."""
    if latmin > latmax or lonmin > lonmax:
        raise ValueError(
            f"Invalid AOI: lat=[{latmin},{latmax}] lon=[{lonmin},{lonmax}]"
        )
    col_min, row_max = latlon_to_tile(latmin, lonmin)  # SW corner -> max row
    col_max, row_min = latlon_to_tile(latmax, lonmax)  # NE corner -> min row
    tiles: list[tuple[int, int]] = []
    for c in range(col_min, col_max + 1):
        for r in range(row_min, row_max + 1):
            tiles.append((c, r))
    return tiles


def tile_filename(col: int, row: int) -> str:
    """`srtm_CC_RR.zip` (zero-padded to 2 digits, matching CGIAR convention)."""
    return f"srtm_{col:02d}_{row:02d}.zip"


def snap_auxdata_dir() -> Path:
    """Resolve SNAP's SRTM 3Sec auxdata folder for the current user.

    Honors `SNAP_USERDIR` env var if set (some sysadmins relocate it),
    otherwise defaults to `~/.snap/auxdata/dem/SRTM 3Sec/`.
    """
    env_userdir = os.environ.get("SNAP_USERDIR")
    if env_userdir:
        base = Path(env_userdir)
    else:
        base = Path.home() / ".snap"
    return base / "auxdata" / "dem" / "SRTM 3Sec"


def download_tile(col: int, row: int, dest: Path,
                  mirrors: list[str] = SRTM_MIRRORS,
                  timeout: int = 60) -> bool:
    """Download a single SRTM tile to `dest`. Tries each mirror in order.
    Returns True on success, False if every mirror fails. Skips download
    if `dest` already exists (caller should pre-check)."""
    name = tile_filename(col, row)
    last_err: Exception | None = None
    for mirror in mirrors:
        url = f"{mirror}/{name}"
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                if response.status != 200:
                    last_err = RuntimeError(
                        f"HTTP {response.status} from {url}")
                    continue
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    while True:
                        chunk = response.read(64 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)
                tmp.rename(dest)
                return True
        except (urllib.error.URLError, OSError) as e:
            last_err = e
            continue
    print(f"WARN: all mirrors failed for {name}: {last_err}",
          file=sys.stderr)
    return False


def parse_project_conf(path: Path) -> dict[str, float]:
    """Read LATMIN / LATMAX / LONMIN / LONMAX from a snap2stamps project.conf.

    Tolerant to whitespace and comments. Format expected (one key per line):
        KEY = VALUE
    """
    keys = {"LATMIN", "LATMAX", "LONMIN", "LONMAX"}
    found: dict[str, float] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.split("#", 1)[0].strip()
        if key in keys:
            try:
                found[key] = float(val)
            except ValueError:
                continue
    missing = keys - found.keys()
    if missing:
        raise SystemExit(
            f"project.conf {path} is missing AOI keys: {sorted(missing)}"
        )
    return found


def _ensure_dir(path: Path) -> None:
    """Create directory tree, robust against Python 3.13 + Windows
    `os.makedirs(exist_ok=True)` regression (WinError 267 on fresh paths)."""
    parts: list[Path] = []
    cur = path
    while not cur.exists():
        parts.append(cur)
        cur = cur.parent
    for p in reversed(parts):
        try:
            os.mkdir(p)
        except FileExistsError:
            pass


def precache(latmin: float, latmax: float, lonmin: float, lonmax: float,
             auxdata: Path | None = None,
             downloader=download_tile) -> dict[str, list[str]]:
    """Ensure all SRTM tiles covering the AOI are present in `auxdata`.

    `downloader` is injectable for testing.
    Returns a report dict: {'present': [...], 'downloaded': [...], 'failed': [...]}.
    """
    if auxdata is None:
        auxdata = snap_auxdata_dir()
    _ensure_dir(auxdata)

    tiles = aoi_tiles(latmin, latmax, lonmin, lonmax)
    report: dict[str, list[str]] = {
        "present": [], "downloaded": [], "failed": [],
    }
    for col, row in tiles:
        name = tile_filename(col, row)
        dest = auxdata / name
        if dest.exists() and dest.stat().st_size > 0:
            report["present"].append(name)
            continue
        ok = downloader(col, row, dest)
        if ok:
            report["downloaded"].append(name)
        else:
            report["failed"].append(name)
    return report


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 1:
        print("Usage: srtm_precache.py <project.conf>", file=sys.stderr)
        return 2
    conf_path = Path(argv[0])
    if not conf_path.is_file():
        print(f"project.conf not found: {conf_path}", file=sys.stderr)
        return 2

    aoi = parse_project_conf(conf_path)
    auxdata = snap_auxdata_dir()
    print(f"AOI: lat=[{aoi['LATMIN']}, {aoi['LATMAX']}], "
          f"lon=[{aoi['LONMIN']}, {aoi['LONMAX']}]")
    print(f"SNAP auxdata: {auxdata}")

    tiles = aoi_tiles(aoi["LATMIN"], aoi["LATMAX"],
                      aoi["LONMIN"], aoi["LONMAX"])
    print(f"SRTM 3Sec tiles required: {len(tiles)}")
    for col, row in tiles:
        print(f"  {tile_filename(col, row)}")

    report = precache(aoi["LATMIN"], aoi["LATMAX"],
                      aoi["LONMIN"], aoi["LONMAX"], auxdata=auxdata)
    print(f"\nDone. present={len(report['present'])}, "
          f"downloaded={len(report['downloaded'])}, "
          f"failed={len(report['failed'])}")
    if report["failed"]:
        print("FAILED tiles (StampsExport may still try to auto-download):")
        for name in report["failed"]:
            print(f"  {name}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
