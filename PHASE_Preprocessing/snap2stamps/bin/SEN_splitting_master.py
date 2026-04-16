### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Author: Jose Manuel Delgado Blasco
# Modified by: Roberto Monti
# Automated Multi-Swath AOI Bounding Box Update

import os
from pathlib import Path
import sys
import glob
import subprocess
import time

inputfile = sys.argv[1]

bar_message='\n#####################################################################\n'

print("## INPUT PARAMETERS ##")
try:
    in_file = open(inputfile, 'r')
    for line in in_file.readlines():
        if "PROJECTFOLDER" in line:
            PROJECT = line.split('=')[1].strip()
            print("\n- PROJECT FOLDER: \n" + PROJECT)
        # We now expect a comma-separated list like SWATHS=IW1,IW2
        if "SWATHS" in line:
            SWATHS_STR = line.split('=')[1].strip()
            swaths = [s.strip() for s in SWATHS_STR.split(',')]
            print("\n- SWATHS TO CHECK: \n" + str(swaths))
        if "POLARISATION" in line:
            POLARISATION = line.split('=')[1].strip()
            print("\n- POLARIZATION: \n" + POLARISATION)
        if "LONMIN" in line:
            LONMIN = line.split('=')[1].strip()
        if "LATMIN" in line:
            LATMIN = line.split('=')[1].strip()
        if "LONMAX" in line:
            LONMAX = line.split('=')[1].strip()
        if "LATMAX" in line:
            LATMAX = line.split('=')[1].strip()
        if "GRAPHSFOLDER" in line:
            GRAPH = line.split('=')[1].strip()
        if "GPTBIN_PATH" in line:
            GPT = line.split('=')[1].strip()
        if "CACHE" in line:
            CACHE = line.split('=')[1].strip()
        if "CPU" in line:
            CPU = line.split('=')[1].strip()
finally:
    in_file.close()

# 1. Define the universal AOI Polygon
polygon = f"POLYGON (({LONMIN} {LATMIN},{LONMAX} {LATMIN},{LONMAX} {LATMAX},{LONMIN} {LATMAX},{LONMIN} {LATMIN}))"
print(f"\n- USING AOI POLYGON: \n{polygon}")

#############################################################################
### TOPSAR Splitting (Assembling) and Apply Orbit section ####
############################################################################
masterfolder=PROJECT+'/master'
splitfolder=PROJECT+'/master'
logfolder=PROJECT+'/logs'
graphfolder=PROJECT+'/graphs'

for f in [splitfolder, logfolder, graphfolder]:
    if not os.path.exists(f):
        os.makedirs(f)

graph2run=PROJECT+'/graphs/SEN_splitgraph2run.xml'
outlog=logfolder+'/SEN_split_proc_stdout.log'
out_file = open(outlog, 'a')
err_file = out_file

print(bar_message)
out_file.write(bar_message)
message='## TOPSAR Splitting and Apply Orbit (Master)\n'
print(message)
out_file.write(message)
print(bar_message)
out_file.write(bar_message)

k = 0
for acdatefolder in sorted(os.listdir(masterfolder)):
    # SAFEGUARD: Skip any hidden files or loose zips, only process directories
    folder_path = os.path.join(masterfolder, acdatefolder)
    if not os.path.isdir(folder_path):
        continue

    k += 1
    print(f'[{k}] Folder: {acdatefolder}')
    out_file.write(f'[{k}] Folder: {acdatefolder}\n')

    files = glob.glob(os.path.join(folder_path, '*.zip'))
    splitmasterfolder = os.path.join(splitfolder, acdatefolder)
    if not os.path.exists(splitmasterfolder):
        os.makedirs(splitmasterfolder)

    # 2. Iterate through each requested swath
    for swath in swaths:
        outputname = f"{acdatefolder}_split_{swath}_Orb.dim"
        output_path = os.path.join(splitmasterfolder, outputname)

        print(f"  -> Testing intersection for {swath}...")

        if len(files) == 1:
            graphxml=GRAPH+'/SEN_master_split_applyorbit.xml'
            with open(graphxml, 'r') as file:
                filedata = file.read()
            filedata = filedata.replace('INPUTFILE', files[0])

        elif len(files) == 2:
            graphxml=GRAPH+'/SEN_master_assemble_split_applyorbit.xml'
            with open(graphxml, 'r') as file:
                filedata = file.read()
            filedata = filedata.replace('INPUTFILE1', files[0])
            filedata = filedata.replace('INPUTFILE2', files[1])

        else:
            print("  -> Invalid number of zip files. Skipping.")
            continue

        # Dynamic variable replacement
        filedata = filedata.replace('IWs', swath)
        filedata = filedata.replace('POLARISATION', POLARISATION)
        filedata = filedata.replace('POLYGON', polygon)
        filedata = filedata.replace('OUTPUTFILE', output_path)

        with open(graph2run, 'w') as file:
            file.write(filedata)

        args = [GPT, graph2run, '-c', CACHE, '-q', CPU]

        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        timeStarted = time.time()
        stdout = process.communicate()[0].decode('utf-8')
        timeDelta = time.time() - timeStarted

        # 3. Graceful error handling for empty intersections
        if process.returncode != 0:
            msg = f"  -> Note: {swath} skipped or failed (Likely no spatial intersection with AOI). Exec time: {timeDelta:.2f}s\n"
            print(msg)
            out_file.write(msg)
        else:
            msg = f"  -> Success: Split master {swath} completed in {timeDelta:.2f}s.\n"
            print(msg)
            out_file.write(msg)

print(bar_message)
out_file.write(bar_message)
out_file.close()
