"""Read a GeoTIFF's declared GDAL NoData value without loading its pixels."""
import math
from pathlib import Path
import struct


def geotiff_nodata(path):
    """Return GDAL_NODATA (TIFF tag 42113), or None when not declared.

    Classic TIFF and BigTIFF, either byte order, are supported. Only the
    first IFD and the short ASCII tag are read; raster size is irrelevant.
    """
    with Path(path).open('rb') as stream:
        def read(count):
            data = stream.read(count)
            if len(data) != count:
                raise ValueError(f'Truncated GeoTIFF header: {path}')
            return data

        marker = read(2)
        if marker not in (b'II', b'MM'):
            raise ValueError(f'External DEM must be a valid GeoTIFF: {path}')
        order = '<' if marker == b'II' else '>'
        version = struct.unpack(order + 'H', read(2))[0]
        if version == 42:
            offset = struct.unpack(order + 'I', read(4))[0]
            count_format, entry_size, inline_size = 'H', 12, 4
        elif version == 43:
            if struct.unpack(order + 'HH', read(4)) != (8, 0):
                raise ValueError('Unsupported BigTIFF offset format.')
            offset = struct.unpack(order + 'Q', read(8))[0]
            count_format, entry_size, inline_size = 'Q', 20, 8
        else:
            raise ValueError(f'Unsupported TIFF version: {version}')
        stream.seek(offset)
        count = struct.unpack(order + count_format, read(struct.calcsize(count_format)))[0]
        if count > 65536:
            raise ValueError('Invalid GeoTIFF tag count.')
        for _ in range(count):
            entry = read(entry_size)
            tag, kind = struct.unpack(order + 'HH', entry[:4])
            if tag != 42113:
                continue
            length = struct.unpack(order + ('I' if version == 42 else 'Q'), entry[4:entry_size-inline_size])[0]
            if kind != 2 or not 1 <= length <= 1024:
                raise ValueError('Invalid GDAL_NODATA GeoTIFF tag.')
            value = entry[-inline_size:]
            if length > inline_size:
                stream.seek(struct.unpack(order + ('I' if version == 42 else 'Q'), value)[0])
                value = read(length)
            raw = value[:length].split(b'\0', 1)[0].decode('ascii').strip()
            result = float(raw)
            if math.isinf(result):
                raise ValueError('Infinite DEM NoData is unsupported; declare a finite value or NaN.')
            return result
    return None


def external_dem_nodata(dem_name, dem_file):
    """Resolve the graph value; retain SNAP's legacy default if undeclared."""
    if dem_name.strip().lower() != 'external dem':
        return '0.0'  # Ignored by SNAP for its built-in DEMs.
    value = geotiff_nodata(dem_file)
    if value is None:
        print('Warning: external DEM has no declared NoData; SNAP default 0.0 is retained.', flush=True)
        return '0.0'
    text = 'NaN' if math.isnan(value) else format(value, '.17g')
    print(f'External DEM NoData read from GeoTIFF: {text}', flush=True)
    return text


def prepare_external_dem(dem_name, dem_file, project):
    """Return a cached SNAP-compatible DEM and its actual NoData value.

    Never edit the supplied DEM or resample elevations. Prefer NoData=0;
    retain valid sea-level pixels by choosing a safe integer sentinel instead
    when zero is a valid elevation. Work in bounded strips, not whole rasters.
    """
    if dem_name.strip().lower() != 'external dem':
        return str(dem_file), '0.0'
    original_nodata = geotiff_nodata(dem_file)
    if original_nodata == 0:
        return str(dem_file), '0.0'
    import hashlib
    import json
    import os
    import tempfile
    try:
        import rasterio
    except ImportError:
        import subprocess
        import sys
        print('Preparing external DEM support in the processing Python interpreter...', flush=True)
        subprocess.run([sys.executable, '-m', 'pip', 'install', '--disable-pip-version-check',
                        '--no-input', '--only-binary=:all:', 'rasterio>=1.4.4,<2'], check=True)
        import importlib
        importlib.invalidate_caches()
        import rasterio
    import numpy as np

    def digest(path):
        h = hashlib.sha256()
        with Path(path).open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(block)
        return h.hexdigest()

    source = Path(dem_file).resolve()
    before = source.stat()
    source_hash = digest(source)
    cache = Path(project) / 'auxiliary' / 'DEM'
    cache.mkdir(parents=True, exist_ok=True)
    output = cache / f'{source.stem}_PHASE_{source_hash[:16]}.tif'
    manifest = output.with_suffix('.json')
    if output.is_file() and manifest.is_file():
        try:
            saved = json.loads(manifest.read_text(encoding='utf-8'))
            if (saved.get('schema') == 1 and saved.get('sourceSHA256') == source_hash
                    and saved.get('outputSHA256') == digest(output)):
                print(f'Using compatible external DEM: {output} (NoData={saved["noData"]:g}).', flush=True)
                return str(output), format(saved['noData'], '.17g')
        except (ValueError, KeyError):
            pass

    fd, temporary_name = tempfile.mkstemp(prefix='.phase-dem-', suffix='.tif', dir=cache)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        with rasterio.open(source) as src:
            if src.count != 1 or src.crs is None or src.crs.to_epsg() != 4326:
                raise ValueError('External DEM must be a single-band WGS84 (EPSG:4326) GeoTIFF; no automatic reprojection is performed.')
            if np.dtype(src.dtypes[0]).kind not in 'iuf':
                raise ValueError('External DEM must contain real numeric elevations.')
            from rasterio.windows import Window

            def windows():
                for row in range(0, src.height, 128):
                    yield Window(0, row, src.width, min(128, src.height-row))

            def read(window):
                values = src.read(1, window=window)
                valid = (src.read_masks(1, window=window) > 0) & np.isfinite(values)
                return values, valid

            candidates = [0, -32768, -9999, -2147483648]
            used = set()
            valid_count = missing_count = 0
            for window in windows():
                values, valid = read(window)
                valid_count += int(valid.sum())
                missing_count += int((~valid).sum())
                for candidate in candidates:
                    if np.any(values[valid] == candidate):
                        used.add(candidate)
            if not valid_count:
                raise ValueError('External DEM has no finite valid elevations.')
            nodata = next((value for value in candidates if value not in used), None)
            if nodata is None:
                raise ValueError('No safe integer NoData sentinel is available without masking valid elevations.')
            if nodata != 0:
                print(f'External DEM has valid zero elevations; preserving them with NoData={nodata}.', flush=True)
            dtype = np.dtype(src.dtypes[0])
            if dtype.kind in 'iu' and not np.iinfo(dtype).min <= nodata <= np.iinfo(dtype).max:
                dtype = np.promote_types(dtype, np.min_scalar_type(nodata))
            profile = src.profile.copy()
            profile.update(driver='GTiff', dtype=dtype.name, nodata=nodata,
                           compress='deflate', BIGTIFF='IF_SAFER')
            # Source layouts can carry a predictor incompatible with promoted data.
            profile.pop('predictor', None)
            with rasterio.open(temporary, 'w', **profile) as dst:
                dst.update_tags(**src.tags())
                dst.update_tags(1, **src.tags(1))
                dst.scales = src.scales
                dst.offsets = src.offsets
                if src.descriptions[0]:
                    dst.set_band_description(1, src.descriptions[0])
                if src.units[0]:
                    dst.set_band_unit(1, src.units[0])
                for window in windows():
                    values, valid = read(window)
                    values = values.astype(dtype, copy=True)
                    values[~valid] = nodata
                    dst.write(values, 1, window=window)
            with rasterio.open(temporary) as dst:
                if (dst.crs != src.crs or dst.transform != src.transform or dst.shape != src.shape
                        or dst.scales != src.scales or dst.offsets != src.offsets):
                    raise RuntimeError('DEM normalization changed raster geometry or elevation scaling.')
                for window in windows():
                    values, valid = read(window)
                    converted = dst.read(1, window=window)
                    if (not np.array_equal(converted[valid], values[valid])
                            or not np.array_equal(dst.read_masks(1, window=window)>0, valid)):
                        raise RuntimeError('DEM normalization changed valid elevations or the missing-data mask.')
        after = source.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError('The supplied DEM changed during normalization. Retry with a stable input file.')
        saved = dict(schema=1, source=str(source), sourceSHA256=source_hash,
                     outputSHA256=digest(temporary), noData=nodata,
                     validPixels=valid_count, missingPixels=missing_count,
                     resamplingApplied=False, validElevationsUnchanged=True)
        os.replace(temporary, output)
        manifest_temp = manifest.with_name(manifest.name+'.tmp')
        manifest_temp.write_text(json.dumps(saved, indent=2)+'\n', encoding='utf-8')
        os.replace(manifest_temp, manifest)
        print(f'Prepared compatible external DEM: {output} (NoData={nodata}; original unchanged).', flush=True)
        return str(output), format(nodata, '.17g')
    finally:
        temporary.unlink(missing_ok=True)
