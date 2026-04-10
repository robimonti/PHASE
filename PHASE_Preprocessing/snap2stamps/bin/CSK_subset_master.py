### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Modified for Flat Architecture and Regex Date Extraction

import os
import sys
import glob
import subprocess
import time
import re

inputfile = sys.argv[1]
bar_message='\n#####################################################################\n'

# Getting configuration variables from inputfile
print("############################# STEP 2 (MASTER) ################################\n")
print("## INPUT PARAMETERS ##")
try:
    with open(inputfile, 'r') as in_file:
        for line in in_file.readlines():
            if "PROJECTFOLDER" in line:
                    PROJECT = line.split('=')[1].strip()
                    print("\n- PROJECT FOLDER: ")
                    print(PROJECT)
            if "GRAPHSFOLDER" in line:
                    GRAPH = line.split('=')[1].strip()
                    print("\n- GRAPH FOLDER: ")
                    print(GRAPH)
            if "GPTBIN_PATH" in line:
                    GPT = line.split('=')[1].strip()
                    print("\n- GPT: ")
                    print(GPT)
            if "LONMIN" in line:
                    LONMIN = line.split('=')[1].strip()
                    print("\n- LONMIN: ")
                    print(LONMIN)
            if "LATMIN" in line:
                    LATMIN = line.split('=')[1].strip()
                    print("\n- LATMIN: ")
                    print(LATMIN)
            if "LONMAX" in line:
                    LONMAX = line.split('=')[1].strip()
                    print("\n- LONMAX: ")
                    print(LONMAX)
            if "LATMAX" in line:
                    LATMAX = line.split('=')[1].strip()
                    print("\n- LATMAX: ")
                    print(LATMAX)
            if "CACHE" in line:
                    CACHE = line.split('=')[1].strip()
            if "CPU" in line:
                    CPU = line.split('=')[1].strip()
except Exception as e:
    print(f"Error reading config: {e}")
    sys.exit(1)

polygon = f'POLYGON (({LONMIN} {LATMIN},{LONMAX} {LATMIN},{LONMAX} {LATMAX},{LONMIN} {LATMAX},{LONMIN} {LATMIN}))'

masterfolder = os.path.join(PROJECT, 'master')
logfolder = os.path.join(PROJECT, 'logs')
graphfolder = os.path.join(PROJECT, 'graphs')

for f in [masterfolder, logfolder, graphfolder]:
    if not os.path.exists(f):
        os.makedirs(f)

graph2run = os.path.join(graphfolder, 'CSK_subsetMgraph2run.xml')
outlog = os.path.join(logfolder, 'CSK_subsetM_proc_stdout.log')

with open(outlog, 'a') as out_file:
    err_file = out_file
    print(bar_message)
    out_file.write(bar_message)
    k = 0
    # Reading from subfolder just like S1
    for acdatefolder in sorted(os.listdir(masterfolder)):
        folder_path = os.path.join(masterfolder, acdatefolder)
        if not os.path.isdir(folder_path): continue

        k += 1
        print(f'[{k}] Folder: {acdatefolder}')
        files = glob.glob(os.path.join(folder_path, '*.h5'))
        if not files: continue

        filename = os.path.basename(files[0])
        # Safe Regex Extraction
        match = re.search(r'_(\d{8})\d{6}_', filename)
        master_date = match.group(1) if match else filename[37:45]

        outputname = f"{master_date}_sub.dim"
        output_path = os.path.join(masterfolder, outputname) # Outputs to /master root like before
        graphxml = os.path.join(GRAPH, 'CSK_AOI_subset.xml')

        with open(graphxml, 'r') as file:
            filedata = file.read()

        filedata = filedata.replace('INPUTFILE', files[0])
        filedata = filedata.replace('OUTPUTFILE', output_path)
        filedata = filedata.replace('POLYGON', polygon)

        with open(graph2run, 'w') as file:
            file.write(filedata)

        args = [GPT, graph2run, '-c', CACHE, '-q', CPU]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        timeStarted = time.time()
        stdout = process.communicate()[0].decode('utf-8', errors='ignore')
        timeDelta = time.time() - timeStarted

        msg = f'[{k}] Finished process in {timeDelta:.2f} seconds.\n'
        print(msg)
        out_file.write(msg)
