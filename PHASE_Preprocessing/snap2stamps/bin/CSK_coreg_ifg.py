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
print("############################# STEP 3 ################################\n")
print("## INPUT PARAMETERS ##")
try:
    with open(inputfile, 'r') as in_file:
        for line in in_file.readlines():
                if "PROJECTFOLDER" in line:
                        PROJECT = line.split('=')[1].strip()
                        print("\n- PROJECT FOLDER: ")
                        print(PROJECT)
                if "DEMNAME" in line:
                        DEMNAME = line.split('=')[1].strip()
                        print("\n- DEM NAME: ")
                        print(DEMNAME)
                if "DEMFILE" in line:
                        DEMFILE = line.split('=')[1].strip()
                        print("\n- DEM FILE: ")
                        print(DEMFILE)
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
                if "NUMGCP" in line:
                        NUMGCP = line.split('=')[1].strip()
                if "CACHE" in line:
                        CACHE = line.split('=')[1].strip()
                if "CPU" in line:
                        CPU = line.split('=')[1].strip()
except Exception as e:
    print(f"Error reading config: {e}")
    sys.exit(1)


subsetfolder = os.path.join(PROJECT, 'subset')
outputcoregfolder = os.path.join(PROJECT, 'coreg')
outputifgfolder = os.path.join(PROJECT, 'ifg')
logfolder = os.path.join(PROJECT, 'logs')
graphsfolder = os.path.join(PROJECT, 'graphs')

for f in [outputcoregfolder, outputifgfolder, logfolder, graphsfolder]:
    if not os.path.exists(f):
        os.makedirs(f)

outlog = os.path.join(logfolder, 'CSK_coreg_ifg_proc_stdout.log')
graphxml = os.path.join(GRAPH, 'CSK_coreg_ifg_computation.xml')
graph2run = os.path.join(graphsfolder, 'CSK_coreg_ifg2run.xml')

with open(outlog, 'a') as out_file:
    err_file = out_file
    print(bar_message)
    out_file.write(bar_message)
    message = '## CSK Coregistration and Interferogram computation started:\n'
    print(message)
    out_file.write(message)
    print(bar_message)
    out_file.write(bar_message)

    # Safely extract Master Date using 8-digit Regex
    master_filename = os.path.basename(MASTER)
    master_match = re.search(r'(\d{8})', master_filename)
    master_date = master_match.group(1) if master_match else "UNKNOWN_MASTER"

    # Flat directory search for slaves (removed the /*/ subfolder requirement)
    slave_files = glob.glob(os.path.join(subsetfolder, '*.dim'))

    for k, dimfile in enumerate(slave_files, start=1):
        slave_filename = os.path.basename(dimfile)
        print(dimfile)
        message = f'[{k}] Processing slave file: {slave_filename}\n'
        print(message)
        out_file.write(message)

        # Safely extract Slave Date using 8-digit Regex
        slave_match = re.search(r'(\d{8})', slave_filename)
        slave_date = slave_match.group(1) if slave_match else "UNKNOWN_SLAVE"

        outputname = f"{master_date}_{slave_date}.dim"

        with open(graphxml, 'r') as file:
            filedata = file.read()

        filedata = filedata.replace('MASTER', MASTER)
        filedata = filedata.replace('SLAVE', dimfile)
        filedata = filedata.replace('DEMNAME', DEMNAME)
        filedata = filedata.replace('DEMFILE', DEMFILE)
        filedata = filedata.replace('NUMGCP', NUMGCP)
        filedata = filedata.replace('OUTPUTCOREGFOLDER', outputcoregfolder)
        filedata = filedata.replace('OUTPUTIFGFOLDER', outputifgfolder)
        filedata = filedata.replace('OUTPUTFILE', outputname)

        with open(graph2run, 'w') as file:
            file.write(filedata)

        args = [GPT, graph2run, '-c', CACHE, '-q', CPU]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        timeStarted = time.time()
        stdout = process.communicate()[0].decode('utf-8', errors='ignore')
        print(f'SNAP STDOUT:\n{stdout}')
        timeDelta = time.time() - timeStarted

        message_finish = f'[{k}] Finished process in {timeDelta:.2f} seconds.\n'
        print(message_finish)
        out_file.write(message_finish)

        if process.returncode != 0:
            err_message = f'Error computing coregistration for subsetted slave {slave_filename}\n'
            err_file.write(err_message)
        else:
            succ_message = f'Coregistration for {slave_filename} successfully completed.\n'
            print(succ_message)
            out_file.write(succ_message)
            print(bar_message)
            out_file.write(bar_message)
