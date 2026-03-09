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

# Getting configuration variables from inputfile
print("############################# STEP 2 ################################\n")
print("## INPUT PARAMETERS ##")
try:
        in_file = open(inputfile, 'r')

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
finally:
        in_file.close()

polygon='POLYGON (('+LONMIN+' '+LATMIN+','+LONMAX+' '+LATMIN+','+LONMAX+' '+LATMAX+','+LONMIN+' '+LATMAX+','+LONMIN+' '+LATMIN+'))'

#############################################################################
########### TOPSAR Splitting (Assembling) and Apply Orbit section ###########
#############################################################################
slavefolder=PROJECT+'/slaves'
subsetfolder=PROJECT+'/subset'
logfolder=PROJECT+'/logs'
graphfolder=PROJECT+'/graphs'
if not os.path.exists(subsetfolder):
                os.makedirs(subsetfolder)
if not os.path.exists(logfolder):
                os.makedirs(logfolder)
if not os.path.exists(graphfolder):
                os.makedirs(graphfolder)

graph2run=PROJECT+'/graphs/CSK_subsetgraph2run.xml'
outlog=logfolder+'/CSK_subset_proc_stdout.log'
out_file = open(outlog, 'a')
err_file=out_file

print(bar_message)
out_file.write(bar_message)
message='## Raster Subset\n'
print(message)
out_file.write(message)
print(bar_message)
out_file.write(bar_message)
k=0
for acdatefolder in sorted(os.listdir(slavefolder)):
    k=k+1
    print('['+str(k)+'] Folder: '+acdatefolder)
    out_file.write('['+str(k)+'] Folder: '+acdatefolder+'\n')
    print(os.path.join(slavefolder, acdatefolder))
    out_file.write(str(os.path.join(slavefolder, acdatefolder))+'\n')
    files = glob.glob(os.path.join(slavefolder, acdatefolder) + '/*.h5')
    print(files)
    out_file.write(str(files)+'\n')
    head, tail = os.path.split(os.path.join(str(files)))
    subslavefolder=subsetfolder+'/'+tail[27:35]
    if not os.path.exists(subslavefolder):
                os.makedirs(subslavefolder)
    outputname=tail[27:35]+'_sub.dim'
    graphxml=GRAPH+'/CSK_AOI_subset.xml'
    # Read in the file
    print('FILE(s) : '+files[0])
    with open(graphxml, 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace('INPUTFILE', files[0])
    filedata = filedata.replace('OUTPUTFILE',subslavefolder+'/'+outputname)
    filedata = filedata.replace('POLYGON',polygon)
    # # Write the file out again
    with open(graph2run, 'w') as file:
        file.write(filedata)
    args = [ GPT, graph2run, '-c', CACHE, '-q', CPU]
    print(args)
    out_file.write(str(args)+'\n')
    # launching the process
    process = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    timeStarted = time.time()
    stdout = process.communicate()[0]
    print('SNAP STDOUT:{}'.format(stdout))
    timeDelta = time.time() - timeStarted                     # Get execution time.
    print(('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.'))
    out_file.write('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.\n')
    if process.returncode != 0 :
        message='Error subsetting image '+str(files)
        err_file.write(message)
    else: 
        message='Subset image '+str(files)+' successfully completed.\n'
        print(message)
        out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)
out_file.close()

