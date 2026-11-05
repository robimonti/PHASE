"""Subset the selected CSK master; report missing inputs or SNAP failures."""
from pathlib import Path
import re
import subprocess
import sys
import time
from xml.sax.saxutils import escape
from phase_subprocess import live_popen, log_console


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
    graph = Path(config['GRAPHSFOLDER']) / 'CSK_AOI_subset.xml'
    polygon = ('POLYGON (({LONMIN} {LATMIN},{LONMAX} {LATMIN},'
               '{LONMAX} {LATMAX},{LONMIN} {LATMAX},{LONMIN} {LATMIN}))').format(**config)
    output = master / f'{date}_sub.dim'
    contents = graph.read_text(encoding='utf-8')
    for key, value in [('INPUTFILE', str(source)), ('OUTPUTFILE', str(output)), ('POLYGON', polygon)]:
        contents = contents.replace(key, escape(value))
    graphs = project / 'graphs'
    graphs.mkdir(parents=True, exist_ok=True)
    generated = graphs / 'CSK_subsetMgraph2run.xml'
    generated.write_text(contents, encoding='utf-8')
    command = [config['GPTBIN_PATH'], str(generated), '-c', config['CACHE'], '-q', config['CPU']]
    print('################ STEP 1: CSK master subsetting ################', flush=True)
    print(f'Master input: {source}\nMaster output: {output}', flush=True)
    print('SNAP command: ' + repr(command), flush=True)
    started = time.time()
    process = live_popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f'SNAP CSK master subsetting failed (exit {process.returncode}). See the SNAP output above.')
    if not output.is_file():
        raise RuntimeError(f'SNAP returned success but did not create the CSK master: {output}')
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
