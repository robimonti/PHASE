"""Subset the selected CSK master; report missing inputs or SNAP failures."""
from pathlib import Path
import re
import subprocess
import sys
import time
from xml.sax.saxutils import escape
from phase_subprocess import log_console
from phase_csk_geometry import subset_product


def subset(config):
    required = ['PROJECTFOLDER', 'GRAPHSFOLDER', 'GPTBIN_PATH', 'LONMIN', 'LONMAX', 'LATMIN', 'LATMAX', 'CACHE', 'CPU']
    missing = [key for key in required if not config.get(key)]
    if missing:
        raise ValueError('Missing CSK master configuration: ' + ', '.join(missing))
    project = Path(config['PROJECTFOLDER'])
    master = project / 'master'
    inputs = sorted(master.rglob('*.h5'))
    if not inputs:
        raise FileNotFoundError(f'No selected CSK master HDF5 file in {master}.')
    dates = {}
    for path in inputs:
        match = re.search(r'_(\d{8})\d{6}_', path.name)
        if not match:
            raise ValueError(f'Cannot determine CSK acquisition date: {path.name}')
        dates.setdefault(match.group(1), []).append(path)
    if len(dates) != 1:
        raise ValueError('Selected master contains multiple acquisition dates: ' + ', '.join(dates))
    date, products = next(iter(dates.items()))
    source = products[0]  # Same first-product policy as the original CSK graph.
    if len(products) > 1:
        print(f'CSK master date {date}: {len(products)} products; subsetting {source.name}.', flush=True)
    output = master / f'{date}_sub.dim'
    graphs = project / 'graphs'
    graphs.mkdir(parents=True, exist_ok=True)
    generated = graphs / 'CSK_subsetMgraph2run.xml'
    print('################ STEP 1: CSK master subsetting ################', flush=True)
    print(f'Master input: {source}\nMaster output: {output}', flush=True)
    started = time.time()
    subset_product(source,output,config,generated)
    print(f'CSK master completed in {time.time() - started:.2f} seconds: {output}', flush=True)


def main():
    config = {}
    for line in Path(sys.argv[1]).read_text(encoding='utf-8-sig').splitlines():
        if line.strip() and not line.lstrip().startswith('#') and '=' in line:
            key, value = line.split('=', 1)
            config[key.strip()] = value.strip()
    if not config.get('PROJECTFOLDER'):
        raise ValueError('Missing PROJECTFOLDER in CSK master configuration.')
    log = Path(config['PROJECTFOLDER']) / 'logs' / 'CSK_subsetM_proc_stdout.log'
    with log_console(log):
        subset(config)


if __name__ == '__main__':
    main()
