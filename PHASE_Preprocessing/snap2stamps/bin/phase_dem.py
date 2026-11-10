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
