### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Modified to prevent recursive `.data` folder scanning

import os
import sys
import glob
import subprocess
import time

inputfile = sys.argv[1]
bar_message = '\n#####################################################################\n'

print("############################# STEP 6 ################################\n")
print("## INPUT PARAMETERS ##")
try:
    with open(inputfile, 'r') as in_file:
        for line in in_file.readlines():
            if "PROJECTFOLDER" in line: PROJECT = line.split('=')[1].strip()
            if "DEMNAME" in line: DEMNAME = line.split('=')[1].strip()
            if "DEMFILE" in line: DEMFILE = line.split('=')[1].strip()
            if "GRAPHSFOLDER" in line: GRAPH = line.split('=')[1].strip()
            if "GPTBIN_PATH" in line: GPT = line.split('=')[1].strip()
            if "TC_COHERENCE" in line: TC_COHERENCE = line.split('=')[1].strip()
            if "CACHE" in line: CACHE = line.split('=')[1].strip()
            if "CPU" in line: CPU = line.split('=')[1].strip()
except Exception as e:
    print(f"Error reading config: {e}")
    sys.exit(1)

coherencefolder = os.path.join(PROJECT, 'coherence')
liafolder = os.path.join(PROJECT, 'lia')
ifgfolder = os.path.join(PROJECT, 'ifg')
logfolder = os.path.join(PROJECT, 'logs')
graphfolder = os.path.join(PROJECT, 'graphs')

for f in [coherencefolder, liafolder, graphfolder]:
    if not os.path.exists(f):
        os.makedirs(f)

outlog = os.path.join(logfolder, 'CSK_tc_proc_stdout.log')
COH_FLAG = int(TC_COHERENCE)
DEMRESAMPLING = 'BILINEAR_INTERPOLATION'

if COH_FLAG == 0:
    graphxml = os.path.join(GRAPH, 'CSK_slave_terrain_correction.xml')
    graph2run = os.path.join(graphfolder, 'CSK_terrain_correction2run.xml')

    with open(outlog, 'a') as out_file:
        err_file = out_file
        print(bar_message)
        out_file.write(bar_message)
        message = '## Terrain correction for coherence and LIA \n'
        print(message)
        out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)

        # Strictly scan top level of ifgfolder using glob to prevent deep dive into .data
        interferogram_files = glob.glob(os.path.join(ifgfolder, '*.dim'))
        
        if not interferogram_files:
            print("No interferogram files found in the ifgfolder.")
            sys.exit(1)

        for k, dimfile in enumerate(interferogram_files, start=1):
            filename = os.path.basename(dimfile)
            print(dimfile)
            message = f'[{k}] Processing interferogram file: {filename}\n'
            print(message)
            out_file.write(message)
            
            tail_ne = filename[:-4]
            folder_path = os.path.join(ifgfolder, f'{tail_ne}.data')
            
            if not os.path.isdir(folder_path):
                continue
                
            coherence_band = None
            for file in os.listdir(folder_path):
                if file.startswith('coh') and file.endswith('.img'):
                    coherence_band = file[:-4]
                    break
                    
            if coherence_band:
                print(f"Coherence band found: {coherence_band}")
            else:
                print(f"No coherence band found for {folder_path}")

            sourceband_lia = 'localIncidenceAngle'
            outputname_coh = os.path.join(coherencefolder, f'{tail_ne}_coh_TC.tif')
            outputname_lia = os.path.join(liafolder, f'{tail_ne}_lia_TC.tif')
            
            with open(graphxml, 'r') as file:
                filedata = file.read()

            filedata = filedata.replace('SLAVE', dimfile)
            filedata = filedata.replace('DEMNAME', DEMNAME)
            filedata = filedata.replace('DEMFILE', DEMFILE)
            filedata = filedata.replace('DEM_RESAMPLING', DEMRESAMPLING)
            filedata = filedata.replace('SOURCEBAND_LIA', sourceband_lia)
            filedata = filedata.replace('OUTPUTFILE_LIA', outputname_lia)
            filedata = filedata.replace('OUTPUTFILE_COH', outputname_coh)
            filedata = filedata.replace('COHBAND', coherence_band)

            with open(graph2run, 'w') as file:
                file.write(filedata)

            args = [GPT, graph2run, '-c', CACHE, '-q', CPU]
            process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            timeStarted = time.time()
            stdout = process.communicate()[0].decode('utf-8', errors='ignore')
            
            print(f'SNAP STDOUT:\n{stdout}')
            timeDelta = time.time() - timeStarted
            msg_finish = f'[{k}] Finished process in {timeDelta:.2f} seconds.\n'
            print(msg_finish)
            out_file.write(msg_finish)

            if process.returncode != 0:
                err_file.write(f'Error terrain correction {filename}\n')
            else:
                succ_msg = f'Terrain correction for coherence and lia, {filename} successfully completed.\n'
                print(succ_msg)
                out_file.write(succ_msg)
                print(bar_message)
                out_file.write(bar_message)
else:
    print(bar_message)
    with open(outlog, 'a') as out_file:
        message = '## Terrain correction for coherence and LIA skipped as TC_COHERENCE is set = 1 \n'
        print(message)
        out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)