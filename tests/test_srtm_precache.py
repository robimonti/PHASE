"""Tests for tools/srtm_precache.py — offline (no network)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import srtm_precache as sp


# ---- latlon_to_tile -----------------------------------------------------

def test_latlon_to_tile_calabria():
    # AOI center used in PHASE test dataset (Calabria)
    col, row = sp.latlon_to_tile(39.12, 16.69)
    assert (col, row) == (40, 5)


def test_latlon_to_tile_origin():
    # Origin: lon=-180, lat=+60 -> tile 1,1
    assert sp.latlon_to_tile(60.0, -180.0) == (1, 1)


def test_latlon_to_tile_meridian():
    # Just east of prime meridian, equator -> tile 37, 12 (lon ~0, lat ~0)
    # row = floor((60-0.5)/5) + 1 = floor(11.9) + 1 = 12
    col, row = sp.latlon_to_tile(0.5, 0.5)
    assert col == 37
    assert row == 12


# ---- aoi_tiles ----------------------------------------------------------

def test_aoi_tiles_calabria_single():
    # AOI ~30 km in Calabria fits in a single 5° tile
    tiles = sp.aoi_tiles(latmin=38.97, latmax=39.27,
                         lonmin=16.51, lonmax=16.86)
    assert tiles == [(40, 5)]


def test_aoi_tiles_spans_two_cols():
    # Crossing a 5° meridian -> 2 cols
    tiles = sp.aoi_tiles(latmin=39.0, latmax=39.5,
                         lonmin=14.5, lonmax=15.5)
    cols = sorted({c for c, _ in tiles})
    assert len(cols) == 2


def test_aoi_tiles_spans_two_rows():
    tiles = sp.aoi_tiles(latmin=34.5, latmax=35.5,
                         lonmin=16.0, lonmax=17.0)
    rows = sorted({r for _, r in tiles})
    assert len(rows) == 2


def test_aoi_tiles_invalid_raises():
    with pytest.raises(ValueError):
        sp.aoi_tiles(latmin=40, latmax=39, lonmin=16, lonmax=17)
    with pytest.raises(ValueError):
        sp.aoi_tiles(latmin=39, latmax=40, lonmin=17, lonmax=16)


# ---- tile_filename ------------------------------------------------------

def test_tile_filename_padding():
    assert sp.tile_filename(1, 5) == "srtm_01_05.zip"
    assert sp.tile_filename(40, 5) == "srtm_40_05.zip"
    assert sp.tile_filename(72, 24) == "srtm_72_24.zip"


# ---- snap_auxdata_dir ---------------------------------------------------

def test_snap_auxdata_dir_default(monkeypatch):
    monkeypatch.delenv("SNAP_USERDIR", raising=False)
    expected = Path.home() / ".snap" / "auxdata" / "dem" / "SRTM 3Sec"
    assert sp.snap_auxdata_dir() == expected


def test_snap_auxdata_dir_env_override(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("SNAP_USERDIR", str(tmp_path))
    assert sp.snap_auxdata_dir() == tmp_path / "auxdata" / "dem" / "SRTM 3Sec"


# ---- parse_project_conf -------------------------------------------------

def test_parse_project_conf_basic(tmp_path: Path):
    f = tmp_path / "project.conf"
    f.write_text(
        "PROJECTFOLDER = /some/path\n"
        "LATMIN = 38.97\n"
        "LATMAX = 39.27\n"
        "LONMIN = 16.51\n"
        "LONMAX = 16.86\n"
        "MASTER = 20210914\n"
    )
    aoi = sp.parse_project_conf(f)
    assert aoi == {"LATMIN": 38.97, "LATMAX": 39.27,
                   "LONMIN": 16.51, "LONMAX": 16.86}


def test_parse_project_conf_with_comments_and_whitespace(tmp_path: Path):
    f = tmp_path / "project.conf"
    f.write_text(
        "# AOI for Calabria test\n"
        "  LATMIN  =   38.97   # south\n"
        "LATMAX=39.27\n"
        "LONMIN = 16.51\n"
        "LONMAX = 16.86 # east\n"
    )
    aoi = sp.parse_project_conf(f)
    assert aoi["LATMIN"] == 38.97
    assert aoi["LONMAX"] == 16.86


def test_parse_project_conf_missing_keys_raises(tmp_path: Path):
    f = tmp_path / "project.conf"
    f.write_text("LATMIN = 1\nLATMAX = 2\n")
    with pytest.raises(SystemExit) as exc:
        sp.parse_project_conf(f)
    assert "LONMIN" in str(exc.value)
    assert "LONMAX" in str(exc.value)


# ---- precache (with injected downloader) --------------------------------

def test_precache_skips_existing_tiles(tmp_path: Path):
    auxdir = tmp_path / "auxdir"
    auxdir.mkdir()
    # Pre-populate with the tile that the AOI requires
    (auxdir / "srtm_40_05.zip").write_bytes(b"fake content")

    calls: list[tuple[int, int]] = []
    def fake_dl(col, row, dest, **kwargs):
        calls.append((col, row))
        return True

    report = sp.precache(38.97, 39.27, 16.51, 16.86,
                        auxdata=auxdir, downloader=fake_dl)
    assert calls == []  # already present, no download
    assert report["present"] == ["srtm_40_05.zip"]
    assert report["downloaded"] == []
    assert report["failed"] == []


def test_precache_downloads_missing(tmp_path: Path):
    auxdir = tmp_path / "auxdir"

    calls: list[tuple[int, int]] = []
    def fake_dl(col, row, dest, **kwargs):
        calls.append((col, row))
        dest.write_bytes(b"downloaded")
        return True

    report = sp.precache(38.97, 39.27, 16.51, 16.86,
                        auxdata=auxdir, downloader=fake_dl)
    assert calls == [(40, 5)]
    assert report["downloaded"] == ["srtm_40_05.zip"]
    assert (auxdir / "srtm_40_05.zip").exists()


def test_precache_reports_failures(tmp_path: Path):
    auxdir = tmp_path / "auxdir"

    def fake_dl(col, row, dest, **kwargs):
        return False  # simulate every mirror failing

    report = sp.precache(38.97, 39.27, 16.51, 16.86,
                        auxdata=auxdir, downloader=fake_dl)
    assert report["failed"] == ["srtm_40_05.zip"]
    assert report["downloaded"] == []


def test_precache_zero_size_treated_as_missing(tmp_path: Path):
    auxdir = tmp_path / "auxdir"
    auxdir.mkdir()
    (auxdir / "srtm_40_05.zip").write_bytes(b"")  # corrupted/empty

    def fake_dl(col, row, dest, **kwargs):
        dest.write_bytes(b"good content")
        return True

    report = sp.precache(38.97, 39.27, 16.51, 16.86,
                        auxdata=auxdir, downloader=fake_dl)
    assert report["downloaded"] == ["srtm_40_05.zip"]


# ---- main CLI -----------------------------------------------------------

def test_main_returns_2_on_missing_arg(capsys):
    rc = sp.main([])
    assert rc == 2
    err = capsys.readouterr().err
    assert "Usage" in err


def test_main_returns_2_on_nonexistent_conf(capsys, tmp_path: Path):
    rc = sp.main([str(tmp_path / "nope.conf")])
    assert rc == 2
    err = capsys.readouterr().err
    assert "not found" in err


def test_main_happy_path(monkeypatch, tmp_path: Path, capsys):
    conf = tmp_path / "project.conf"
    conf.write_text(
        "LATMIN = 38.97\nLATMAX = 39.27\n"
        "LONMIN = 16.51\nLONMAX = 16.86\n"
    )
    auxdir = tmp_path / "auxdir"
    monkeypatch.setattr(sp, "snap_auxdata_dir", lambda: auxdir)

    def fake_dl(col, row, dest, **kwargs):
        dest.write_bytes(b"x")
        return True
    monkeypatch.setattr(sp, "download_tile", fake_dl)

    rc = sp.main([str(conf)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "AOI:" in out
    assert "srtm_40_05.zip" in out


# ---- Live network test (opt-in) -----------------------------------------

@pytest.mark.skipif(
    "PHASE_RUN_NETWORK_TESTS" not in __import__("os").environ,
    reason="set PHASE_RUN_NETWORK_TESTS=1 to enable live SRTM download test"
)
def test_live_download_real_srtm_tile(tmp_path: Path):
    """Real download from CGIAR. Validates the URL/mirror logic against
    the actual server. Marked opt-in to keep CI fast and offline-safe."""
    dest = tmp_path / "srtm_40_05.zip"
    ok = sp.download_tile(40, 5, dest, timeout=120)
    assert ok, "Could not download srtm_40_05.zip from any mirror"
    assert dest.stat().st_size > 1000, "Downloaded file too small to be valid"
