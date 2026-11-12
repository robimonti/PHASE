"""GeoTIFF NoData propagation in the actual CSK graph-generation scripts."""
import importlib.util
from pathlib import Path
import runpy
import struct
import sys
import xml.etree.ElementTree as ET
import pytest

BIN = Path(__file__).resolve().parents[1] / 'PHASE_Preprocessing/snap2stamps/bin'
spec = importlib.util.spec_from_file_location('phase_dem', BIN / 'phase_dem.py')
dem = importlib.util.module_from_spec(spec); spec.loader.exec_module(dem)


def tiff(path, value, order='<', big=False):
    marker = b'II' if order == '<' else b'MM'
    raw = (value + '\0').encode() if value is not None else b''
    if big:
        header = marker + struct.pack(order+'HHHQ', 43, 8, 0, 16)
        entry = struct.pack(order+'HHQ', 42113, 2, len(raw))
        entry += raw.ljust(8, b'\0') if len(raw) <= 8 else struct.pack(order+'Q', 52)
        body = struct.pack(order+'Q', int(value is not None)) + (entry if raw else b'') + struct.pack(order+'Q', 0)
    else:
        header = marker + struct.pack(order+'HI', 42, 8)
        entry = struct.pack(order+'HHI', 42113, 2, len(raw))
        entry += raw.ljust(4, b'\0') if len(raw) <= 4 else struct.pack(order+'I', 26)
        body = struct.pack(order+'H', int(value is not None)) + (entry if raw else b'') + struct.pack(order+'I', 0)
    path.write_bytes(header + body + (raw if len(raw) > (8 if big else 4) else b''))


@pytest.mark.parametrize('order', ['<', '>'])
@pytest.mark.parametrize('big', [False, True])
@pytest.mark.parametrize('value', ['-3.4028234663852886e+38', '-9999', '0'])
def test_declared_nodata_is_read_without_pixel_data(tmp_path, order, big, value):
    path = tmp_path / 'dem.tif'; tiff(path, value, order, big)
    assert dem.geotiff_nodata(path) == float(value)
    assert float(dem.external_dem_nodata('External DEM', path)) == float(value)


def test_undeclared_value_retains_documented_legacy_default(tmp_path, capsys):
    path = tmp_path / 'dem.tif'; tiff(path, None)
    assert dem.geotiff_nodata(path) is None
    assert dem.external_dem_nodata('External DEM', path) == '0.0'
    assert 'no declared NoData' in capsys.readouterr().out
    assert dem.external_dem_nodata('SRTM 1Sec HGT', 'unused') == '0.0'


def test_corrupt_dem_is_not_silently_used(tmp_path):
    path = tmp_path / 'dem.tif'; path.write_bytes(b'II*\0')
    with pytest.raises(ValueError, match='Truncated'):
        dem.external_dem_nodata('External DEM', path)


@pytest.mark.parametrize('exit_code', [0, 7])
def test_actual_coreg_script_uses_declared_nodata_in_both_operators(tmp_path, monkeypatch, phase_root, exit_code):
    for folder in ['master', 'subset', 'graphs']:(tmp_path / folder).mkdir()
    master = tmp_path / 'master/20251011_sub.dim'; master.touch()
    (tmp_path / 'subset/20250504_sub.dim').touch()
    sys.path.insert(0,str(BIN))
    from phase_csk_geometry import save_products
    save_products(tmp_path/'subset',[tmp_path/'subset/20250504_sub.dim'])
    path = tmp_path / 'dem & elevations.tif'; real_dem(path, [[500, 600], [-3.4028234663852886e+38, 800]])
    config = {'PROJECTFOLDER':tmp_path, 'DEMNAME':'External DEM', 'DEMFILE':path, 'MASTER':master,
              'GRAPHSFOLDER':phase_root/'PHASE_Preprocessing/snap2stamps/graphs', 'GPTBIN_PATH':'unused',
              'NUMGCP':1000, 'CACHE':'128M', 'CPU':1}
    conf = tmp_path / 'project.conf';conf.write_text('\n'.join(f'{k}={v}' for k,v in config.items()))
    monkeypatch.syspath_prepend(str(BIN))
    import phase_subprocess
    class Process:
        returncode=exit_code
        def communicate(self):
            if not exit_code:
                for folder in ['coreg','ifg']:(tmp_path/folder/'20251011_20250504.dim').touch()
            return (b'controlled GPT completed', None)
    monkeypatch.setattr(phase_subprocess, 'live_popen', lambda *a,**kw:Process())
    monkeypatch.setattr(sys, 'argv', ['CSK_coreg_ifg.py',str(conf)])
    if exit_code:
        with pytest.raises(RuntimeError, match='exit 7'):
            runpy.run_path(str(BIN/'CSK_coreg_ifg.py'),run_name='__main__')
        assert 'controlled GPT completed' in (tmp_path/'logs/CSK_coreg_ifg_proc_stdout.log').read_text()
    else:
        runpy.run_path(str(BIN/'CSK_coreg_ifg.py'),run_name='__main__')
    graph = ET.parse(tmp_path/'graphs/CSK_coreg_ifg2run.xml')
    values=[float(n.text) for n in graph.findall('.//externalDEMNoDataValue')]
    assert values == [0.0]*2
    prepared = graph.find('.//externalDEMFile').text
    assert Path(prepared) != path and Path(prepared).is_file()


def real_dem(path, values, nodata=-3.4028234663852886e38, dtype='float32', crs='EPSG:4326'):
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
    data = np.asarray(values, dtype=dtype)
    with rasterio.open(path, 'w', driver='GTiff', count=1, width=data.shape[1], height=data.shape[0],
                       dtype=dtype, crs=crs, transform=from_origin(8.7, 46.1, .0001, .0001), nodata=nodata) as r:
        r.write(data, 1)
        r.update_tags(AREA_OR_POINT='Area', example='retain')
        r.scales = (2.,)
        r.offsets = (3.,)
    return data


@pytest.mark.parametrize('missing', [-3.4028234663852886e38, -9999., float('nan')])
def test_normalization_preserves_geometry_valid_pixels_and_cache(tmp_path, missing):
    import numpy as np
    import rasterio
    source = tmp_path/'source.tif'
    original = real_dem(source, [[500, 600], [missing, 800]], nodata=missing)
    before = source.read_bytes()
    output, nodata = dem.prepare_external_dem('External DEM', source, tmp_path)
    assert float(nodata) == 0
    with rasterio.open(source) as a, rasterio.open(output) as b:
        valid = a.read_masks(1)>0
        assert np.array_equal(b.read(1)[valid], original[valid])
        assert np.array_equal(b.read_masks(1)>0, valid)
        assert (b.transform, b.crs, b.shape, b.scales, b.offsets) == (a.transform,a.crs,a.shape,a.scales,a.offsets)
        assert b.tags()['example'] == 'retain'
    assert source.read_bytes() == before
    timestamp = Path(output).stat().st_mtime_ns
    assert dem.prepare_external_dem('External DEM', source, tmp_path)[0] == output
    assert Path(output).stat().st_mtime_ns == timestamp
    # A damaged cached file must be regenerated rather than trusted.
    Path(output).write_bytes(b'damaged cache')
    assert dem.prepare_external_dem('External DEM', source, tmp_path)[0] == output
    with rasterio.open(output) as b:assert b.nodata == 0


def test_valid_zero_elevations_are_never_masked(tmp_path):
    import numpy as np
    import rasterio
    source = tmp_path/'source.tif'
    real_dem(source, [[0, 100], [-9999, 800]], nodata=-9999)
    output, nodata = dem.prepare_external_dem('External DEM', source, tmp_path)
    assert float(nodata) == -32768
    with rasterio.open(output) as b:
        assert b.read(1)[0,0] == 0 and b.read_masks(1)[0,0] == 255
        assert b.read(1)[1,0] == -32768 and b.read_masks(1)[1,0] == 0


def test_declared_zero_and_builtin_dems_need_no_conversion(tmp_path):
    source=tmp_path/'source.tif';real_dem(source, [[0, 100]], nodata=0)
    assert dem.prepare_external_dem('External DEM', source, tmp_path) == (str(source),'0.0')
    assert dem.prepare_external_dem('SRTM 1Sec HGT', '', tmp_path) == ('','0.0')
    assert not (tmp_path/'auxiliary').exists()


def test_unsupported_crs_is_not_silently_reprojected(tmp_path):
    source=tmp_path/'source.tif';real_dem(source, [[1,2]], nodata=-9999, crs='EPSG:32632')
    with pytest.raises(ValueError, match='WGS84'):
        dem.prepare_external_dem('External DEM', source, tmp_path)
    assert not list(tmp_path.rglob('.phase-dem-*.tif'))
