import os
from pathlib import Path
import sys
import subprocess
import shlex
import time

inputfile = sys.argv[1]

bar_message = '\n#####################################################################\n'

# Getting configuration variables from inputfile
print("############################# STEP 6 ################################\n")
print("## INPUT PARAMETERS ##")
try:
    in_file = open(inputfile, 'r')

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
        if "GRAPHSFOLDER" in line:
            GRAPH = line.split('=')[1].strip()
            print("\n- GRAPH FOLDER: ")
            print(GRAPH)
        if "GPTBIN_PATH" in line:
            GPT = line.split('=')[1].strip()
            print("\n- GPT: ")
            print(GPT)
        if "TC_COHERENCE" in line:
            TC_COHERENCE = line.split('=')[1].strip()
            print("\n- COHERENCE FLAG: ")
            print(TC_COHERENCE)
        if "CACHE" in line:
            CACHE = line.split('=')[1].strip()
        if "CPU" in line:
            CPU = line.split('=')[1].strip()
finally:
    in_file.close()

#############################################################################
########### Coherence and LIA terrain correction and export (tif) ###########
#############################################################################
coherencefolder = PROJECT + '/coherence'
liafolder = PROJECT + '/lia'
ifgfolder = PROJECT + '/ifg'
logfolder = PROJECT + '/logs'
graphfolder = PROJECT + '/graphs'
if not os.path.exists(coherencefolder):
    os.makedirs(coherencefolder)
if not os.path.exists(liafolder):
    os.makedirs(liafolder)
if not os.path.exists(graphfolder):
    os.makedirs(graphfolder)

outlog = logfolder + '/CSK_tc_proc_stdout.log'

COH_FLAG = int(TC_COHERENCE)
DEMRESAMPLING = 'BILINEAR_INTERPOLATION'

if COH_FLAG == 0:
    graphxml = GRAPH + '/CSK_slave_terrain_correction.xml'
    print(graphxml)
    graph2run = PROJECT + '/graphs/CSK_terrain_correction2run.xml'

    print(bar_message)
    with open(outlog, 'a') as out_file:
        err_file = out_file

        out_file.write(bar_message)
        message = '## Terrain correction for coherence and LIA \n'
        print(message)
        out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)

        # List interferogram files in the ifgfolder
        interferogram_files = []
        for root, dirs, files in os.walk(ifgfolder):
            for file in files:
                if file.endswith('.dim'):
                    interferogram_files.append(os.path.join(root, file))

        print("Interferogram files found:")
        for file in interferogram_files:
            print(file)

        # Check if the ifgfolder contains any files
        if not interferogram_files:
            print("No interferogram files found in the ifgfolder.")
            sys.exit(1)

        k = 0
        for dimfile in interferogram_files:
            print(dimfile)
            k = k + 1
            head, tail = os.path.split(dimfile)
            message = '[' + str(k) + '] Processing interferogram file :' + tail + '\n'
            print(message)
            out_file.write(message)
            
            tail_ne = tail[:-4]
            
            # Iterate through folders in 'ifg' directory
            folder_name = tail_ne + '.data'
            folder_path = os.path.join(ifgfolder, folder_name)
            if not os.path.isdir(folder_path):
                continue  # Skip if not a directory
            
            if folder_name.endswith('.data'):
                coherence_band = None
                for file in os.listdir(folder_path):
                    if file.startswith('coh') and file.endswith('.img'):
                        coherence_band = file[:-4]
                        break
            
                if coherence_band:
                    print(("Coherence band found:", coherence_band))
                else:
                    print(("No coherence band found for", folder_name))

            sourceband_lia = 'localIncidenceAngle'
            outputname_coh = tail_ne + '_coh_TC.tif'
            outputname_lia = tail_ne + '_lia_TC.tif'
            with open(graphxml, 'r') as file:
                filedata = file.read()

            # Replace the target string
            filedata = filedata.replace('SLAVE', dimfile)
            filedata = filedata.replace('DEMNAME', DEMNAME)
            filedata = filedata.replace('DEMFILE', DEMFILE)
            filedata = filedata.replace('DEM_RESAMPLING', DEMRESAMPLING)
            filedata = filedata.replace('SOURCEBAND_LIA', sourceband_lia)
            filedata = filedata.replace('OUTPUTFILE_LIA', outputname_lia)
            filedata = filedata.replace('OUTPUTFILE_COH', outputname_coh)
            filedata = filedata.replace('OUTPUTLIAFOLDER', liafolder)
            filedata = filedata.replace('OUTPUTCOHFOLDER', coherencefolder)
            filedata = filedata.replace('COHBAND', coherence_band)

            # Write the file out again
            with open(graph2run, 'w') as file:
                file.write(filedata)

            args = [GPT, graph2run, '-c', CACHE, '-q', CPU]

            # Launch the processing
            process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            timeStarted = time.time()
            stdout = process.communicate()[0]
            print(('SNAP STDOUT:{}'.format(stdout)))
            timeDelta = time.time() - timeStarted  # Get execution time.
            print(('[' + str(k) + '] Finished process in ' + str(timeDelta) + ' seconds.'))
            out_file.write('[' + str(k) + '] Finished process in ' + str(timeDelta) + ' seconds.\n')

            if process.returncode != 0:
                message = 'Error terrain correction ' + str(tail)
                err_file.write(message)
            else:
                message = 'Terrain correction for coherence and lia, ' + str(tail) + ' successfully completed.\n'
                print(message)
                out_file.write(message)
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

