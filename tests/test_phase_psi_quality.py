"""Scientific invariants for read-only preliminary PSI screening."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

BIN = Path(__file__).resolve().parents[1]/'PHASE_Preprocessing/snap2stamps/bin'
sys.path.insert(0, str(BIN))
spec = importlib.util.spec_from_file_location('phase_psi_quality', BIN/'phase_psi_quality.py')
quality = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quality)


def test_amplitude_dispersion_is_invariant_to_acquisition_gain_and_masks():
    grid = quality.grid_definition([9, 9.001, 45, 45.001], 'CSK')
    lon = np.linspace(9, 9.001, 400); lat = np.linspace(45, 45.001, 400)
    amplitude = np.linspace(1, 10, 400)
    slcs = [(amplitude*g*np.exp(1j*np.arange(400))).astype('complex64') for g in (1, 2, 4, 3, 7, 5)]
    data = quality.stack_statistics(slcs, lon, lat, grid, progress=lambda *a, **k: None)
    observed = data['samples'] > 0
    np.testing.assert_allclose(data['coverage'][observed], 1)
    np.testing.assert_allclose(data['stability'][observed], 1)
    assert data['samples'].sum() == 400
    slcs[0][:100] = 0
    data = quality.stack_statistics(slcs, lon, lat, grid, progress=lambda *a, **k: None)
    assert np.nanmin(data['coverage']) < 1
    assert np.all(data['amplitudeSamples'] <= data['samples'])


def test_unknown_geometry_never_creates_favourable_class_and_risk_dominates():
    data = dict(samples=np.full(4, 100), coverage=np.ones(4), stability=np.ones(4))
    classes, reasons = quality.classify(data, np.array([np.nan, 0, .3, 0]), np.array([.8, .8, .8, np.nan]), 10)
    np.testing.assert_array_equal(classes, [0, 1, 3, 0])
    assert reasons[0] & 32 and reasons[2] & 8
    classes, _ = quality.classify(data, np.zeros(4), np.ones(4), 2)
    assert not classes.any()


def test_back_facing_terrain_uses_surface_normal_and_satellite_direction():
    lon, lat = np.meshgrid(np.linspace(0, .001, 5), np.linspace(0, .001, 5))
    q = np.zeros((25, 10));q[:, 0] = lon.ravel();q[:, 1] = lat.ravel()
    q[:, 7] = 7000000  # Overhead satellite at the equator: flat terrain faces it.
    back, valid = quality.back_facing_terrain(q, 5, 5)
    assert valid.all() and not back.any()
    q[:, 7] = 5000000  # Opposite direction: same surface normal faces away.
    back, valid = quality.back_facing_terrain(q, 5, 5)
    assert valid.all() and back.all()


def test_geometry_helper_is_portable_and_matches_its_source():
    import hashlib
    folder = BIN.parent/'java';name = 'PhaseQualityGeometry'
    manifest = json.loads((folder/(name+'.build.json')).read_text())
    source = (folder/(name+'.java')).read_bytes().replace(b'\r\n', b'\n')
    compiled = (folder/(name+'.class')).read_bytes()
    assert hashlib.sha256(source).hexdigest() == manifest['sourceSHA256']
    assert hashlib.sha256(compiled).hexdigest() == manifest['classSHA256']
    assert int.from_bytes(compiled[6:8], 'big') == 52


def test_shared_map_controls_do_not_diverge_between_modules():
    ui = BIN.parents[1]
    assert (ui/'phase_preprocessing_beta_ui/quality.js').read_bytes() == (ui/'phase_stamps_beta_ui/quality.js').read_bytes()


def test_cell_aggregation_does_not_wrap_invalid_coordinates_or_reverse_north():
    grid = quality.grid_definition([9, 9.001, 45, 45.001], 'CSK')
    cells, valid = quality.cell_indices(np.array([9, 9.001, np.nan, 0, 181]), np.array([45.001, 45, 45, 0, 45]), grid)
    assert cells[0] == 0
    assert cells[1] == grid['rows']*grid['columns']-1
    np.testing.assert_array_equal(valid, [True, True, False, False, False])
    with pytest.raises(ValueError): quality.grid_definition([179, -179, 45, 46], 'CSK')


def fixture(tmp_path):
    export = tmp_path/'INSAR_20250101'; (export/'rslc').mkdir(parents=True); (export/'geo').mkdir()
    width, height = 40, 30
    (export/'rslc/20250101.rslc.par').write_text(f'range_samples: {width}\nazimuth_lines: {height}\n')
    x, y = np.meshgrid(np.linspace(9, 9.001, width), np.linspace(45.001, 45, height))
    x.astype('>f4').tofile(export/'geo/20250101.lon'); y.astype('>f4').tofile(export/'geo/20250101.lat')
    for day, gain in zip(range(1, 7), (1, 2, 3, 5, 4, 6)):
        np.full((height, width), gain+gain*1j, dtype='>c8').tofile(export/f'rslc/2025010{day}.rslc')
    folder = tmp_path/'ifg/20250101_20250102.data'; folder.mkdir(parents=True)
    folder.with_suffix('.dim').write_text('fixture')
    (folder/'coh.hdr').write_text(f'samples = {width}\nlines = {height}\ndata type = 4\nbyte order = 1\n')
    np.full((height, width), .8, dtype='>f4').tofile(folder/'coh.img')
    return dict(project=str(tmp_path), master='20250101', bounds=[9, 9.001, 45, 45.001], sensor='CSK')


def test_full_export_is_read_only_georeferenced_and_keeps_unavailable_evidence(tmp_path):
    rasterio = pytest.importorskip('rasterio')
    request = fixture(tmp_path)
    inputs = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    report = quality.run(request)
    assert report['acquisitionCount'] == 6 and report['interferogramCount'] == 1
    assert report['geometry']['status'] == 'unavailable'
    assert all(value == inputs[path] for path, value in ((p, p.read_bytes()) for p in inputs))
    folder = tmp_path/'INSAR_20250101/quality'
    assert json.loads((folder/'quality.json').read_text())['algorithm'] == quality.ALGORITHM
    with rasterio.open(folder/'coverage.tif') as ds:
        assert ds.crs.to_epsg() == 4326
        assert ds.transform.e < 0
        np.testing.assert_allclose(ds.read(1)[np.isfinite(ds.read(1))], 1)
    with rasterio.open(folder/'geometry.tif') as ds: assert np.isnan(ds.read(1)).all()
    with rasterio.open(folder/'suitability.tif') as ds: assert not np.any(ds.read(1) == 1)


def test_long_stack_keeps_file_mappings_bounded_and_blocks_independent(tmp_path):
    paths = []
    for k in range(5):
        path = tmp_path/f'{k}.rslc';np.full(20, k+1j, dtype='>c8').tofile(path);paths.append(path)
    with quality.BoundedStack(limit=2) as cache:
        saved = cache.raster(paths[0])[:]
        for path in paths:
            samples = cache.raster(path)[:]
            assert len(cache.maps) <= 2
            assert samples.size == 20
        np.testing.assert_array_equal(saved, np.full(20, 1j))
    assert not cache.maps


def test_concurrent_input_change_never_publishes_a_current_report(tmp_path, monkeypatch):
    pytest.importorskip('rasterio');request = fixture(tmp_path)
    original = quality.stack_statistics
    def changed(*args, **kwargs):
        data = original(*args, **kwargs)
        with (tmp_path/'INSAR_20250101/geo/20250101.lon').open('ab') as stream: stream.write(b'x')
        return data
    monkeypatch.setattr(quality, 'stack_statistics', changed)
    with pytest.raises(RuntimeError, match='changed during'):
        quality.run(request)
    assert not (tmp_path/'INSAR_20250101/quality/quality.json').exists()


def test_incompatible_orbit_grid_cannot_produce_geometry_values(tmp_path, monkeypatch):
    rasterio = pytest.importorskip('rasterio');request = fixture(tmp_path)
    import phase_csk_geometry
    from rasterio.transform import from_bounds
    dem = tmp_path/'dem.tif'
    with rasterio.open(dem, 'w', driver='GTiff', width=2, height=2, count=1,
                       dtype='float32', crs='EPSG:4326', nodata=0,
                       transform=from_bounds(9, 45, 9.001, 45.001, 2, 2)) as dst:
        dst.write(np.full((2, 2), 500, dtype='float32'), 1)
    image = tmp_path/'ifg/20250101_20250102.data/coh.img'
    np.full((30, 40), 500, dtype='>f4').tofile(image.parent/'elevation.img')
    request.update(demName='External DEM', demFile=str(dem), gpt='unused')
    monkeypatch.setattr(phase_csk_geometry, 'snap_runtime', lambda _: ('unused', []))
    monkeypatch.setattr(phase_csk_geometry, 'classpath_jar', lambda *a: tmp_path/'unused.jar')
    def bad_model(command, **kwargs):
        points = np.loadtxt(command[-2], delimiter=',', ndmin=2)
        np.savetxt(command[-1], np.full((len(points), 2), 1e6), delimiter=',')
    monkeypatch.setattr(quality.subprocess, 'run', bad_model)
    lon = np.memmap(tmp_path/'INSAR_20250101/geo/20250101.lon', dtype='>f4', mode='r')
    lat = np.memmap(tmp_path/'INSAR_20250101/geo/20250101.lat', dtype='>f4', mode='r')
    folder = tmp_path/'screen';folder.mkdir()
    values, info, _ = quality.geometry_statistics(request, lon, lat, 40, 30,
        quality.grid_definition(request['bounds'], 'CSK'), [(image, '>f4', 1)], folder)
    assert np.isnan(values).all() and info['status'] == 'unavailable'
    assert 'does not reproduce' in info['message']
