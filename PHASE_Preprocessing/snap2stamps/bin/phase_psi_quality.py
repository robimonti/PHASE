"""Preliminary PSI suitability, computed after SNAP export without changing it.

All raster statistics are streamed. Geographic screening cells are descriptive,
not PS probabilities. Geometry is used only after its orbit model reproduces
the actual exported radar grid; unsupported/missing evidence stays unknown.
"""
import argparse
from collections import OrderedDict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

ALGORITHM = 'psi-screening-v1'
BLOCK = 262144
THRESHOLDS = dict(minAcquisitions=5, minRadarSamples=20, amplitudeDispersion=0.4,
                  weakStableFraction=0.01, strongStableFraction=0.05,
                  weakCoverage=0.8, strongCoverage=0.95,
                  weakCoherence=0.2, strongCoherence=0.4,
                  moderateGeometryRisk=0.05, highGeometryRisk=0.2)
LAYERS = {
    'suitability': ('Preliminary PSI suitability', 'Screening classes, not PS probabilities.'),
    'geometry': ('Terrain mapping risk', 'Fraction of projected terrain samples with mapping ambiguity at the analysis scale; not a validated layover/shadow mask.'),
    'stability': ('Amplitude stability', 'Fraction of sufficiently observed radar pixels with normalised amplitude dispersion below 0.4.'),
    'coverage': ('Stack coverage', 'Mean fraction of acquisitions with finite, nonzero complex samples.'),
    'coherence': ('Interferogram coherence', 'Median of spatial cell means across the retained interferograms; not StaMPS temporal coherence.'),
}


class BoundedStack:
    """Limit open mappings as well as RAM for long acquisition stacks."""
    def __init__(self, limit=16):
        self.limit = limit
        self.maps = OrderedDict()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        for raster in self.maps.values():
            raster._mmap.close()
        self.maps.clear()

    def raster(self, path):
        owner = self
        class Raster:
            def __getitem__(self, selection):
                if path not in owner.maps:
                    if len(owner.maps) >= owner.limit:
                        _, old = owner.maps.popitem(last=False); old._mmap.close()
                    owner.maps[path] = np.memmap(path, dtype='>c8', mode='r')
                owner.maps.move_to_end(path)
                # Own the returned block: evicting a mapping cannot invalidate
                # a caller's samples while another acquisition is being read.
                return np.array(owner.maps[path][selection], copy=True)
        return Raster()


def fingerprint(path):
    p = Path(path).resolve()
    s = p.stat()
    return dict(path=str(p), bytes=s.st_size, modified=s.st_mtime)


def grid_definition(bounds, sensor):
    west, east, south, north = map(float, bounds)
    if not (-180 <= west < east <= 180 and -85 < south < north < 85 and east-west < 180):
        raise ValueError('Quality screening requires a non-dateline AOI between 85 S and 85 N.')
    metres = np.array([(east-west)*111320*math.cos(math.radians((south+north)/2)),
                       (north-south)*111320])
    # Respect the spatial support of the existing coherence windows. The map
    # does not claim pixel-scale PS prediction or DEM-scale output accuracy.
    minimum = 100 if str(sensor).upper() == 'SEN' else 50
    cell = max(minimum, float(max(metres))/96)
    cols, rows = np.clip(np.ceil(metres/cell).astype(int), 1, 96)
    return dict(bounds=list(bounds), rows=int(rows), columns=int(cols),
                cellSizeMeters=[float(metres[0]/cols), float(metres[1]/rows)],
                order='Rows north to south; columns west to east')


def cell_indices(lon, lat, grid):
    w, e, s, n = grid['bounds']; nx, ny = grid['columns'], grid['rows']
    valid = np.isfinite(lon) & np.isfinite(lat) & ~((lon == 0) & (lat == 0))
    valid &= (lon >= w) & (lon <= e) & (lat >= s) & (lat <= n)
    x = np.floor(np.clip(np.nan_to_num((lon-w)/(e-w)), 0, 1)*nx).astype(np.int64)
    y = np.floor(np.clip(np.nan_to_num((n-lat)/(n-s)), 0, 1)*ny).astype(np.int64)
    return np.minimum(y, ny-1)*nx + np.minimum(x, nx-1), valid


def aggregate(values, cells, valid, size):
    return np.bincount(cells[valid], weights=np.asarray(values)[valid], minlength=size)


def ratio(numerator, denominator):
    out = np.full(np.shape(numerator), np.nan, dtype=float)
    np.divide(numerator, denominator, out=out, where=np.asarray(denominator) > 0)
    return out


def read_grid(export, master):
    export = Path(export)
    pars = sorted(export.glob(f'*slc/{master}.*slc.par'))
    if not pars:
        available = sorted(p.name for p in export.parent.glob('INSAR_*') if p.is_dir())
        extra = ' Available export folders: '+', '.join(available)+'.' if available else ''
        raise ValueError(f'Master {master}: radar parameters are missing in {export}. '
                         f'Expected rslc/{master}.rslc.par or slc/{master}.slc.par.'+extra)
    text = pars[0].read_text(errors='replace')
    def number(key):
        match = re.search(r'^\s*'+key+r'\s*:?\s*(\d+)', text, re.M)
        if not match:
            raise ValueError('Missing '+key+' in master radar parameters.')
        return int(match[1])
    width, height = number('range_samples'), number('azimuth_lines')
    paths = []
    for suffix in ('lon', 'lat'):
        exact = export/'geo'/f'{master}.{suffix}'
        candidates = sorted(export.glob('geo/*.'+suffix))
        if exact.is_file(): paths.append(exact)
        elif len(candidates) == 1: paths.append(candidates[0])
        elif not candidates: paths.append(None)
        else: raise ValueError(f'Ambiguous {suffix} geolocation rasters in {export}/geo for master {master}.')
    if any(p is None or p.stat().st_size != 4*width*height for p in paths):
        raise ValueError('Export geolocation sizes do not match the radar grid.')
    return width, height, paths, pars[0]


def stack_statistics(slcs, lon, lat, grid, progress=print):
    size = grid['rows']*grid['columns']; total = len(lon)
    counts = np.zeros(size); coverage = np.zeros(size); usable = np.zeros(size); stable = np.zeros(size)
    scales = []
    for k, slc in enumerate(slcs):
        summed, count = 0., 0
        for start in range(0, total, BLOCK):
            stop = min(total, start+BLOCK); _, geo = cell_indices(lon[start:stop], lat[start:stop], grid)
            a = np.abs(slc[start:stop]).astype(float); good = geo & np.isfinite(a) & (a > 0)
            summed += a[good].sum(); count += int(good.sum())
        scales.append(summed/count if count else np.nan)
        progress(f'Quality screening: normalising acquisition {k+1}/{len(slcs)}', flush=True)
    minimum = max(THRESHOLDS['minAcquisitions'], math.ceil(.8*len(slcs)))
    for start in range(0, total, BLOCK):
        stop = min(total, start+BLOCK); cells, geo = cell_indices(lon[start:stop], lat[start:stop], grid)
        n = np.zeros(stop-start); mean = n.copy(); m2 = n.copy()
        for slc, scale in zip(slcs, scales):
            if not np.isfinite(scale) or scale <= 0:
                continue
            a = np.abs(slc[start:stop]).astype(float)/scale
            good = np.isfinite(a) & (a > 0)
            n[good] += 1; delta = a[good]-mean[good]; mean[good] += delta/n[good]
            m2[good] += delta*(a[good]-mean[good])
        enough = geo & (n >= minimum) & (mean > 0)
        # Population standard deviation, as in the StaMPS amplitude screening.
        dispersion = ratio(np.sqrt(ratio(m2, n)), mean)
        counts += np.bincount(cells[geo], minlength=size)
        coverage += aggregate(n/max(1, len(slcs)), cells, geo, size)
        usable += np.bincount(cells[enough], minlength=size)
        good = enough & (dispersion < THRESHOLDS['amplitudeDispersion'])
        stable += np.bincount(cells[good], minlength=size)
    return dict(samples=counts, coverage=ratio(coverage, counts), stability=ratio(stable, usable),
                amplitudeSamples=usable, acquisitionMeanAmplitudes=scales)


def coherence_products(project, master, dates, width, height):
    """Only current-stack interferograms on the exported radar grid qualify."""
    folder = Path(project)/'ifg'
    inventory = folder/'phase_products.json'
    if inventory.is_file():
        from phase_csk_geometry import product_paths
        headers = [Path(p) for p in product_paths(folder)]
    else:
        headers = sorted(folder.glob('*.dim'))
    out = []
    for dim in headers:
        pair = re.findall(r'\d{8}', dim.stem)
        if len(pair) != 2 or master not in pair or not set(pair).issubset(dates):
            continue
        for hdr in sorted(dim.with_suffix('.data').glob('coh*.hdr')):
            text = hdr.read_text(errors='replace')
            fields = {m[1].strip().lower(): m[2].strip() for m in re.finditer(r'^\s*([^=\r\n]+)\s*=\s*([^\r\n]+)', text, re.M)}
            if (fields.get('samples') != str(width) or fields.get('lines') != str(height)
                    or fields.get('data type') != '4' or fields.get('byte order') not in ('0', '1')):
                continue
            image = hdr.with_suffix('.img')
            if image.is_file() and image.stat().st_size == 4*width*height:
                days = abs((datetime.strptime(pair[0], '%Y%m%d')-datetime.strptime(pair[1], '%Y%m%d')).days)
                out.append((image, '>f4' if fields['byte order'] == '1' else '<f4', days))
                break
    return out


def coherence_statistics(products, lon, lat, grid, pair_stack=None):
    size = grid['rows']*grid['columns']; values = []
    for path, dtype, _ in products:
        raster = np.memmap(path, dtype=dtype, mode='r'); summed = np.zeros(size); counts = np.zeros(size)
        for start in range(0, len(lon), BLOCK):
            stop = min(len(lon), start+BLOCK); cells, valid = cell_indices(lon[start:stop], lat[start:stop], grid)
            a = raster[start:stop]; valid &= np.isfinite(a) & (a >= 0) & (a <= 1)
            if pair_stack is not None:
                pair = re.findall(r'\d{8}', path.parent.stem)
                if len(pair) != 2 or not all(date in pair_stack for date in pair):
                    continue
                for date in pair:
                    complex_samples = pair_stack[date][start:stop]
                    valid &= np.isfinite(complex_samples) & (np.abs(complex_samples) > 0)
            summed += aggregate(a, cells, valid, size); counts += np.bincount(cells[valid], minlength=size)
        values.append(ratio(summed, counts))
        del raster
    if not values:
        return np.full(size, np.nan)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmedian(values, axis=0)


def bilinear(raster, x, y, width, height):
    good = np.isfinite(x) & np.isfinite(y) & (x >= 0) & (x < width-1) & (y >= 0) & (y < height-1)
    xx = np.clip(np.nan_to_num(x), 0, width-2); yy = np.clip(np.nan_to_num(y), 0, height-2)
    i, j = np.floor(xx).astype(int), np.floor(yy).astype(int); fx, fy = xx-i, yy-j
    v = np.asarray(raster).reshape(height, width)
    out = (v[j, i]*(1-fx)*(1-fy)+v[j, i+1]*fx*(1-fy)+v[j+1, i]*(1-fx)*fy+v[j+1, i+1]*fx*fy)
    return np.where(good, out, np.nan)


def back_facing_terrain(points, rows, cols):
    """Local DEM normal versus satellite LOS; not a cast-shadow ray tracer."""
    lon, lat, h = np.deg2rad(points[:, 0]), np.deg2rad(points[:, 1]), points[:, 2]
    e2 = 6.6943799901413165e-3
    radius = 6378137/np.sqrt(1-e2*np.sin(lat)**2)
    ground = np.c_[(radius+h)*np.cos(lat)*np.cos(lon),
                   (radius+h)*np.cos(lat)*np.sin(lon),
                   (radius*(1-e2)+h)*np.sin(lat)]
    surface = ground.reshape(rows, cols, 3)
    north, east = np.gradient(surface, axis=(0, 1))
    normal = np.cross(east, north).reshape(-1, 3)
    dot = np.einsum('ij,ij->i', normal, points[:, 7:10]-ground)
    return dot <= 0, np.isfinite(dot) & (np.linalg.norm(normal, axis=1) > 0)


def geometry_statistics(request, lon, lat, width, height, grid, products, folder):
    """Report ambiguity only if SNAP's model matches the actual export grid."""
    size = grid['rows']*grid['columns']; unknown = np.full(size, np.nan)
    info = dict(status='unavailable', message='No compatible terrain/orbit model was available.')
    sources = []
    try:
        from phase_csk_geometry import snap_runtime, classpath_jar
        from phase_dem import prepare_external_dem
        candidates = []
        # The elevation band must be on precisely the exported grid. It supplies
        # ellipsoidal heights for the orbit/grid validation, not ground truth.
        for image, _, _ in products:
            data = image.parent; elevation = data/'elevation.img'
            if elevation.is_file() and elevation.stat().st_size == 4*width*height:
                candidates.append((data.with_suffix('.dim'), elevation))
        if not candidates:
            return unknown, info, sources
        source, elevation = candidates[0]
        dem, nodata = prepare_external_dem(request['demName'], request.get('demFile', ''), request['project'])
        java, jars = snap_runtime(request['gpt']); helper = Path(__file__).resolve().parent.parent/'java'
        name = 'PhaseQualityGeometry'; build = json.loads((helper/(name+'.build.json')).read_text())
        for suffix, key in (('.class', 'classSHA256'), ('.java', 'sourceSHA256')):
            raw = (helper/(name+suffix)).read_bytes()
            if suffix == '.java': raw = raw.replace(b'\r\n', b'\n')
            if hashlib.sha256(raw).hexdigest() != build[key]:
                raise ValueError('Quality geometry helper checksum mismatch.')
        z = np.memmap(elevation, dtype='>f4', mode='r')
        ids = np.unique(np.linspace(0, len(lon)-1, min(512, len(lon))).astype(int))
        valid = np.isfinite(lon[ids]) & np.isfinite(lat[ids]) & np.isfinite(z[ids]) & (np.abs(z[ids]) < 12000)
        valid &= (np.abs(lat[ids]) < 85) & (np.abs(lon[ids]) <= 180) & ~((lon[ids] == 0) & (lat[ids] == 0))
        ids = ids[valid]
        if len(ids) < 12:
            raise ValueError('Too few exported samples to validate the orbit/grid model.')
        validation = folder/'geometry_validation.csv'; projected = folder/'terrain_projection.csv'; check = folder/'geometry_check.csv'
        np.savetxt(validation, np.c_[lon[ids], lat[ids], z[ids]], delimiter=',', fmt='%.12f')
        del z
        # Bounded analysis grid. This is screening at a declared scale, not an
        # exhaustive DEM-resolution visibility/occlusion simulation.
        gx = min(513, max(33, grid['columns']*8+1)); gy = min(513, max(33, grid['rows']*8+1))
        command = [java, '-Xmx1G', '-cp', str(classpath_jar(helper, jars, request['project'])), name,
                   str(source), request['demName'], dem, nodata, *map(str, grid['bounds']),
                   str(gx), str(gy), str(projected), str(validation), str(check)]
        with (folder/'geometry.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
        xy = np.loadtxt(check, delimiter=',', ndmin=2)
        error = np.linalg.norm(xy-np.c_[ids % width, ids//width], axis=1)
        if not np.all(np.isfinite(error)) or np.quantile(error, .95) > 1 or np.max(error) > 3:
            raise ValueError('The terrain/orbit model does not reproduce the exported radar grid; geometry remains unassessed.')
        q = np.loadtxt(projected, delimiter=',', skiprows=1)
        x, y, fx, fy = [q[:, k].reshape(gy, gx) for k in (3, 4, 5, 6)]
        def determinant(a, b):
            ay, ax = np.gradient(a); by, bx = np.gradient(b)
            return ax*by-ay*bx
        det, reference = determinant(x, y), determinant(fx, fy)
        folds = (det*reference <= 0).ravel()
        back_facing, normal_valid = back_facing_terrain(q, gy, gx)
        olon = bilinear(lon, q[:, 3], q[:, 4], width, height)
        olat = bilinear(lat, q[:, 3], q[:, 4], width, height)
        scale = 111320*math.cos(math.radians(sum(grid['bounds'][2:])/2))
        distance = np.hypot((olon-q[:, 0])*scale, (olat-q[:, 1])*111320)
        spacing = max((grid['bounds'][1]-grid['bounds'][0])*scale/(gx-1),
                      (grid['bounds'][3]-grid['bounds'][2])*111320/(gy-1))
        limit = max(20, 2*spacing)
        cells, within = cell_indices(q[:, 0], q[:, 1], grid)
        supported = within & np.isfinite(distance) & np.isfinite(det.ravel()) & (np.abs(reference.ravel()) > 1e-9) & normal_valid
        denominator = np.bincount(cells[within], minlength=size)
        count = np.bincount(cells[supported], minlength=size)
        risk = folds | back_facing | (distance > limit)
        numerator = np.bincount(cells[supported & risk], minlength=size)
        values = ratio(numerator, count); values[ratio(count, denominator) < .8] = np.nan
        info = dict(status='screened', message='Terrain folding, back-facing slopes and coordinate inconsistency. Cast shadows beyond the analysis window are not assessed.',
                    screeningSpacingMeters=spacing, mismatchThresholdMeters=limit,
                    validation95PercentilePixels=float(np.quantile(error, .95)),
                    limitation='DEM-based screening, not a validated layover/shadow mask or an independent geolocation accuracy assessment.')
        sources = [source, elevation]
        if Path(dem).is_file(): sources.append(Path(dem))
        if Path(request.get('demFile', '')).is_file(): sources.append(Path(request['demFile']))
        return values, info, sources
    except Exception as exc:
        info['message'] = str(exc) if isinstance(exc, ValueError) else 'Terrain screening unavailable. See geometry.log for details.'
        info['diagnostic'] = str(exc)
        return unknown, info, sources


def classify(data, geometry, coherence, acquisitions):
    t = THRESHOLDS; n = len(geometry); classes = np.zeros(n, dtype=np.uint8); reasons = np.zeros(n, dtype=np.uint8)
    enough = (data['samples'] >= t['minRadarSamples']) & (acquisitions >= t['minAcquisitions'])
    reasons[~enough] |= 16
    reasons[data['coverage'] < t['weakCoverage']] |= 1
    reasons[data['stability'] < t['weakStableFraction']] |= 2
    reasons[coherence < t['weakCoherence']] |= 4
    reasons[geometry >= t['highGeometryRisk']] |= 8
    available = np.isfinite(geometry) & np.isfinite(coherence) & np.isfinite(data['stability'])
    reasons[~available] |= 32
    classes[enough & available] = 2
    strong = ((data['coverage'] >= t['strongCoverage']) & (data['stability'] >= t['strongStableFraction'])
              & (coherence >= t['strongCoherence']) & (geometry < t['moderateGeometryRisk']))
    classes[enough & available & strong] = 1
    # Poor observed evidence can still be reported if another indicator is
    # missing; missing evidence can never manufacture a favourable class.
    classes[enough & ((reasons & 15) != 0)] = 3
    return classes, reasons


def json_matrix(array, grid):
    return [[None if not np.isfinite(v) else float(v) for v in row]
            for row in np.asarray(array).reshape(grid['rows'], grid['columns'])]


def export_layers(folder, grid, layers):
    import rasterio
    from rasterio.transform import from_bounds
    w, e, s, n = grid['bounds']
    for name, values in layers.items():
        path = folder/(name+'.tif')
        with rasterio.open(path, 'w', driver='GTiff', width=grid['columns'], height=grid['rows'],
                           count=1, dtype='float32', crs='EPSG:4326',
                           transform=from_bounds(w, s, e, n, grid['columns'], grid['rows']),
                           nodata=float('nan'), compress='deflate') as dst:
            dst.write(np.asarray(values, dtype='float32').reshape(grid['rows'], grid['columns']), 1)
            dst.update_tags(algorithm=ALGORITHM, description=LAYERS[name][1],
                            interpretation='Preliminary screening only; does not change StaMPS selection.')


def run(request):
    project, master = Path(request['project']), str(request['master'])
    export = project/('INSAR_'+master); output = export/'quality'
    width, height, paths, par = read_grid(export, master); grid = grid_definition(request['bounds'], request.get('sensor', ''))
    files = sorted(export.glob('rslc/*.rslc'))
    if not files:
        files = sorted(export.glob('slc/*.slc'))
    if not files or any(f.stat().st_size != 8*width*height for f in files):
        raise ValueError('Quality screening requires complete complex float32 SNAP SLC exports on one radar grid.')
    if master not in {p.stem for p in files} or any(not re.fullmatch(r'\d{8}', p.stem) for p in files):
        raise ValueError('The quality stack must contain its master and date-named SNAP acquisitions.')
    output.mkdir(parents=True, exist_ok=True)
    lon, lat = [np.memmap(p, dtype='>f4', mode='r') for p in paths]
    with BoundedStack() as stack:
        slcs = [stack.raster(p) for p in files]
        dates = {f.stem for f in files}; products = coherence_products(project, master, dates, width, height)
        initial = [fingerprint(p) for p in [par, *paths, *files, *(f for f, _, _ in products)]]
        # Stage a complete report; never leave a half-written manifest labelled current.
        manifest = output/'quality.json'; manifest.unlink(missing_ok=True)
        with tempfile.TemporaryDirectory(prefix='phase_quality_', dir=output) as temporary:
            temp = Path(temporary)
            data = stack_statistics(slcs, lon, lat, grid)
            coherence = coherence_statistics(products, lon, lat, grid, dict(zip((f.stem for f in files), slcs)))
            print('Quality screening: checking terrain geometry', flush=True)
            geometry, info, geo_sources = geometry_statistics(request, lon, lat, width, height, grid, products, temp)
            classes, reasons = classify(data, geometry, coherence, len(files))
            layers = dict(suitability=classes.astype(float), geometry=geometry, stability=data['stability'], coverage=data['coverage'], coherence=coherence)
            layers['suitability'][data['samples'] == 0] = np.nan
            export_layers(temp, grid, layers)
            if any(fingerprint(p['path']) != p for p in initial):
                raise RuntimeError('Preprocessing inputs changed during quality screening. Rebuild the map after processing finishes.')
            sources = [par, *paths, *files, *(f for f, _, _ in products), *geo_sources]
            report = dict(schema=1, algorithm=ALGORITHM, generatedAt=datetime.now(timezone.utc).isoformat(),
                          grid=grid, master=master, acquisitionCount=len(files), interferogramCount=len(products),
                          masterResolution=request.get('masterResolution', {}),
                          stackFiles=[p.relative_to(export).as_posix() for p in files],
                          interferogramBaselinesDays=[days for _, _, days in products],
                          thresholds=THRESHOLDS, geometry=info,
                          interpretation='Preliminary screening classes, not PS probabilities. No changes to StaMPS selection. Coherence depends on spatial support and temporal/perpendicular baselines.',
                          amplitudeNormalisation='Each acquisition divided by its mean finite nonzero amplitude inside the AOI; population dispersion across available dates.',
                          acquisitionNormalisation=[dict(date=f.stem, meanAmplitude=float(scale) if np.isfinite(scale) else None)
                                                    for f, scale in zip(files, data['acquisitionMeanAmplitudes'])],
                          classes=['Not evaluated', 'Favourable', 'Intermediate', 'Unfavourable'],
                          reasonBits=[dict(bit=bit, message=message) for bit, message in (
                              (1, 'Incomplete stack coverage'), (2, 'Few amplitude-stable samples'),
                              (4, 'Low interferogram coherence'), (8, 'Terrain mapping risk'),
                              (16, 'Insufficient observations'), (32, 'Indicator unavailable'))],
                          sources=[dict(fingerprint(p),
                                        relativePath=os.path.relpath(p, project),
                                        requiredForFreshness=p in [par, *paths, *files])
                                   for p in dict.fromkeys(sources)],
                          sampleCounts=json_matrix(data['samples'], grid), reasonCodes=json_matrix(reasons, grid),
                          layers={name: dict(title=LAYERS[name][0], description=LAYERS[name][1], values=json_matrix(a, grid)) for name, a in layers.items()})
            for f in temp.glob('*.tif'): os.replace(f, output/f.name)
            if (temp/'geometry.log').is_file(): os.replace(temp/'geometry.log', output/'geometry.log')
            staged = output/'quality.json.tmp'; staged.write_text(json.dumps(report, allow_nan=False), encoding='utf-8'); os.replace(staged, manifest)
        print('PSI quality screening saved: '+str(manifest), flush=True)
        return report



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--request', required=True)
    args = parser.parse_args()
    try:
        run(json.loads(Path(args.request).read_text(encoding='utf-8')))
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(1, f'Quality screening unavailable: {exc}\n')
