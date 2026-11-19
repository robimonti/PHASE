"""Prepare and verify GeoTIFF support in the interpreter running PHASE."""
import argparse
import importlib
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

RASTERIO_PACKAGE = 'rasterio>=1.3.10,<2'


def configure_bundled_data():
    # MATLAB/Conda may export PROJ_LIB/GDAL_DATA for another GDAL build.
    # Wheels carry matching databases: use them only inside this process.
    spec = importlib.util.find_spec('rasterio')
    if spec is None or spec.origin is None:
        return
    package = Path(spec.origin).parent
    if (package / 'proj_data/proj.db').is_file():
        for name in ('PROJ_DATA', 'PROJ_LIB'):
            os.environ[name] = str(package / 'proj_data')
    if (package / 'gdal_data').is_dir():
        os.environ['GDAL_DATA'] = str(package / 'gdal_data')


def check_geotiff():
    """Exercise GDAL, CRS handling and the writer, rather than just importing."""
    import numpy as np
    from rasterio.io import MemoryFile
    from rasterio.transform import from_bounds
    data = np.array([[1, 2], [3, 4]], dtype='float32')
    with MemoryFile() as memory:
        with memory.open(driver='GTiff', width=2, height=2, count=1,
                         dtype='float32', crs='EPSG:4326', nodata=float('nan'),
                         transform=from_bounds(9, 45, 9.01, 45.01, 2, 2),
                         compress='deflate') as dst:
            dst.write(data, 1)
        with memory.open() as src:
            if src.crs.to_epsg() != 4326 or not np.array_equal(src.read(1), data):
                raise RuntimeError('GeoTIFF round-trip validation failed.')


def ensure_geotiff(install=True):
    configure_bundled_data()
    missing = False
    try:
        importlib.import_module('rasterio')
        importlib.import_module('numpy')
    except ModuleNotFoundError as exc:
        if exc.name not in {'rasterio', 'numpy'}:
            raise RuntimeError(f'Broken GeoTIFF dependency in {sys.executable}: {exc}') from exc
        missing = True
    except (ImportError, OSError) as exc:
        raise RuntimeError(f'Cannot load GeoTIFF support in {sys.executable}: {exc}') from exc
    if missing:
        if not install:
            raise RuntimeError(f'GeoTIFF support is missing in {sys.executable}.')
        print(f'Preparing GeoTIFF support (Rasterio/NumPy) in {sys.executable}...', flush=True)
        try:
            subprocess.run([sys.executable, '-m', 'pip', 'install',
                            '--disable-pip-version-check', '--no-input',
                            '--only-binary=:all:', RASTERIO_PACKAGE], check=True)
            # Verify a fresh process so imports made before pip cannot leave
            # stale native-library state in the quality-screening process.
            subprocess.run([sys.executable, str(Path(__file__).resolve())], check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            raise RuntimeError(f'Could not prepare GeoTIFF support in {sys.executable}. '
                               'Check the pip output, internet connection and write permissions; '
                               'then retry Build quality map.') from exc
        importlib.invalidate_caches()
        configure_bundled_data()
        return
    try:
        check_geotiff()
    except Exception as exc:
        raise RuntimeError(f'GeoTIFF support is installed but failed its write/read check '
                           f'in {sys.executable}: {exc}') from exc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ensure', action='store_true')
    args = parser.parse_args()
    try:
        ensure_geotiff(install=args.ensure)
    except RuntimeError as exc:
        parser.exit(1, f'{exc}\n')
    print(f'GeoTIFF support ready: Rasterio/NumPy write and read verified ({sys.executable}).')
