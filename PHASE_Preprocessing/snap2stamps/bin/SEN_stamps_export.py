### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Author: Jose Manuel Delgado Blasco
# Modified by: Roberto Monti
# Date: 13/10/2018
# Version: 2.0

# Step 1 : preparing slaves in folder structure
# Step 2 : TOPSAR Splitting (Assembling) and Apply Orbit
# Step 3 : Coregistration and Interferogram generation
# Step 4 : StaMPS export

# Added option for CACHE and CPU specification by user
# Planned support for DEM selection and ORBIT type selection


import os
from pathlib import Path
import sys
import glob
import subprocess
import shlex
import time
inputfile = sys.argv[1]
bar_message='\n#####################################################################\n'

# SNAP 13 fix: StampsExport has SRTM 3Sec hardcoded for the lat/lon
# geocoding files, regardless of the DEM the user picked for coreg/ifg.
# When gpt cannot auto-download (offline, mirror down, or — on SNAP 13 with
# Sentinel-1C/1D inputs — silent failure), StampsExport produces partially
# corrupted geo files and StaMPS hangs mid-PSI without a clear error.
# Pre-cache the required SRTM tiles before running gpt to eliminate the
# auto-download path entirely.
_phase_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_phase_root / 'tools'))
try:
    import srtm_precache as _srtm
    _aoi = _srtm.parse_project_conf(Path(inputfile))
    print(bar_message)
    print('## SRTM 3Sec pre-cache (SNAP StampsExport DEM workaround) ##')
    _report = _srtm.precache(
        _aoi['LATMIN'], _aoi['LATMAX'], _aoi['LONMIN'], _aoi['LONMAX']
    )
    print(f"  present={len(_report['present'])}, "
          f"downloaded={len(_report['downloaded'])}, "
          f"failed={len(_report['failed'])}")
    if _report['failed']:
        print(f"WARN: {len(_report['failed'])} SRTM tile(s) could not be "
              f"pre-cached: {_report['failed']}")
        print('      StampsExport may still try to auto-download these. '
              'If StaMPS hangs mid-PSI, manually fetch them into '
              "%USERPROFILE%/.snap/auxdata/dem/SRTM 3Sec/ before re-running.")
except SystemExit as _e:
    # AOI keys missing in project.conf: not fatal here; user may run
    # StampsExport with a different DEM that doesn't trigger SRTM auto-download.
    print(f"WARN: SRTM pre-cache skipped — {_e}")
except Exception as _e:
    print(f"WARN: SRTM pre-cache skipped due to unexpected error: {_e}")

# SNAP 13 NPE fix: detect .dim products written by an older SNAP than the
# gpt that's about to read them. The TPG serialization changed between
# SNAP 9 and SNAP 13 — feeding legacy .dim into SNAP 13's StampsExport
# triggers a silent NullPointerException and StaMPS later hangs mid-PSI.
# Warn the user up-front to re-run the pipeline from .SAFE.zip with one
# consistent SNAP version.
try:
    import snap_dim_version_check as _sdv

    _conf_text = Path(inputfile).read_text(encoding="utf-8")
    _coreg_dir = None
    _ifg_dir = None
    _gpt = None
    _project = None
    for _ln in _conf_text.splitlines():
        if "=" not in _ln:
            continue
        _k, _, _v = _ln.partition("=")
        _k = _k.strip(); _v = _v.split("#", 1)[0].strip()
        if _k == "PROJECTFOLDER":
            _project = Path(_v)
        elif _k == "GPTBIN_PATH":
            _gpt = Path(_v)
    if _project is not None and _gpt is not None:
        _coreg_dir = _project / "coreg"
        _ifg_dir = _project / "ifg"
        _dims = []
        for _d in (_coreg_dir, _ifg_dir):
            if _d.is_dir():
                _dims.extend(sorted(_d.glob("*.dim")))
        if _dims:
            _ver_report = _sdv.check_compatibility(_dims, _gpt)
            _sdv.warn_if_mismatch(_ver_report)
except Exception as _e:
    print(f"WARN: SNAP version-mismatch check skipped: {_e}")

# Getting configuration variables from inputfile
print("############################# STEP 4 ################################\n")
print("## INPUT PARAMETERS ##")
try:
        in_file = open(inputfile, 'r')
        for line in in_file.readlines():
                if "PROJECTFOLDER" in line:
                        PROJECT = line.split('=')[1].strip()
                        print("\n- PROJECT FOLDER: ")
                        print(PROJECT)
                if "MASTER" in line:
                        MASTER = line.split('=')[1].strip()
                        print("\n- MASTER: ")
                        print(MASTER)
                if "GRAPHSFOLDER" in line:
                        GRAPH = line.split('=')[1].strip()
                        print("\n- GRAPH FOLDER: ")
                        print(GRAPH)
                if "GPTBIN_PATH" in line:
                        GPT = line.split('=')[1].strip()
                        print("\n- GPT: ")
                        print(GPT)
                if "CACHE" in line:
                        CACHE = line.split('=')[1].strip()
                if "CPU" in line:
                        CPU = line.split('=')[1].strip()
finally:
        in_file.close()

###################################################################################
################################ StaMPS PSI export ################################
###################################################################################
coregfolder=PROJECT+'/coreg'
ifgfolder=PROJECT+'/ifg'
head, tail = os.path.split(MASTER)
outputexportfolder=PROJECT+'/INSAR_'+tail[17:25]
logfolder=PROJECT+'/logs'

if not os.path.exists(outputexportfolder):
                os.makedirs(outputexportfolder)
if not os.path.exists(logfolder):
                os.makedirs(logfolder)

outlog=logfolder+'/SEN_export_proc_stdout.log'
out_file = open(outlog, 'a')
err_file=out_file
graphxml=GRAPH+'/SEN_export.xml'
graph2run=PROJECT+'/graphs/SEN_export2run.xml'
print(bar_message)
out_file.write(bar_message)
message='## StaMPS PSI export started:\n'
print(message)
out_file.write(message)
print(bar_message)
out_file.write(bar_message)
k=0
for dimfile in glob.iglob(coregfolder + '/*.dim'):
    head, tail = os.path.split(os.path.join(coregfolder, dimfile))
    k=k+1
    message='['+str(k)+'] Exporting pair: master-slave pair '+tail+'\n'
    ifgdim = Path(ifgfolder+'/'+tail)
    print(ifgdim)
    if ifgdim.is_file():
        print(message)
        out_file.write(message)
        with open(graphxml, 'r') as file :
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('COREGFILE',dimfile)
        filedata = filedata.replace('IFGFILE', str(ifgdim))
        filedata = filedata.replace('OUTPUTFOLDER',outputexportfolder)
        # Write the file out again
        with open(graph2run, 'w') as file:
            file.write(filedata)
        args = [ GPT, graph2run, '-c', CACHE, '-q', CPU]
        print(args)
        # Launching process
        process = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        timeStarted = time.time()
        stdout = process.communicate()[0]
        print('SNAP STDOUT:{}'.format(stdout))
        timeDelta = time.time() - timeStarted                     # Get execution time.
        print(('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.'))
        out_file.write('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.\n')
        if process.returncode != 0 :
           message='Error exporting '+str(tail)+'\n'
           err_file.write(message)
        else:
           message='Stamps export of '+str(tail)+' successfully completed.\n'
           print(message)
           out_file.write(message)
           print(bar_message)
           out_file.write(bar_message)
out_file.close()
err_file.close()
