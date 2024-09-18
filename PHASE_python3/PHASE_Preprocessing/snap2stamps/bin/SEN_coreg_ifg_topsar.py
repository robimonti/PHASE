### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Author Jose Manuel Delgado Blasco
# Modified by: Roberto Monti
# Date: 11/05/2023
# Version: 2.1

# Step 1 : preparing slaves in folder structure
# Step 2 : TOPSAR Splitting (Assembling) and Apply Orbit
# Step 3 : Coregistration and Interferogram generation
# Step 4 : StaMPS export

# Added option for CACHE and CPU specification by user
# Added support for DEM selection


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
print("############################# STEP 3 ################################\n")
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
                if "COREGDTMNAME" in line:
                        COREGDTMNAME = line.split('=')[1].strip()
                        print("\n- DEM NAME COREG: ")
                        print(COREGDTMNAME)
                if "COREGDTMFILE" in line:
                        COREGDTMFILE = line.split('=')[1].strip()
                        print("\n- DEM FILE COREG: ")
                        print(COREGDTMFILE)
                if "DEMRESAMPLING" in line:
                        DEMRESAMPLING = line.split('=')[1].strip()
                        print("\n- DEM RESAMPLING METHOD: ")
                        print(DEMRESAMPLING)
                if "IW1" in line:
                        IW = line.split('=')[1].strip()
                        print("\n- SWATH: ")
                        print(IW)
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
                if "FIRSTBURSTINDEX" in line:
                    FIRSTBURSTINDEX = line.split('=')[1].strip()
                    print("\n- FIRST BURST: ")
                    print(FIRSTBURSTINDEX)
                if "LASTBURSTINDEX" in line:
                    LASTBURSTINDEX = line.split('=')[1].strip()
                    print("\n- LAST BURST: ")
                    print(LASTBURSTINDEX)
                if "FSTBURSTMST" in line:
                    FSTBURSTMST = line.split('=')[1].strip()
                if "LSTBURSTMST" in line:
                    LSTBURSTMST = line.split('=')[1].strip()
                if "LONMIN" in line:
                    LONMIN = line.split('=')[1].strip()
                if "LATMIN" in line:
                        LATMIN = line.split('=')[1].strip()
                if "LONMAX" in line:
                        LONMAX = line.split('=')[1].strip()
                if "LATMAX" in line:
                        LATMAX = line.split('=')[1].strip()
                if "CACHE" in line:
                    CACHE = line.split('=')[1].strip()
                if "CPU" in line:
                    CPU = line.split('=')[1].strip()
finally:
        in_file.close()

polygon='POLYGON (('+LONMIN+' '+LATMIN+','+LONMAX+' '+LATMIN+','+LONMAX+' '+LATMAX+','+LONMIN+' '+LATMAX+','+LONMIN+' '+LATMIN+'))'

FIRSTBURST = int(FSTBURSTMST)
LASTBURST = int(LSTBURSTMST)
Nburst = LASTBURST - FIRSTBURST

######################################################################################
################### TOPSAR Coregistration and Interferogram formation ################
######################################################################################
slavesplittedfolder=PROJECT+'/split'
outputcoregfolder=PROJECT+'/coreg'
outputifgfolder=PROJECT+'/ifg'
logfolder=PROJECT+'/logs'
graphsfolder=PROJECT+'/graphs'
if not os.path.exists(outputcoregfolder):
                os.makedirs(outputcoregfolder)
if not os.path.exists(outputifgfolder):
                os.makedirs(outputifgfolder)
if not os.path.exists(logfolder):
                os.makedirs(logfolder)
if not os.path.exists(graphsfolder):
                os.makedirs(graphsfolder)
                
outlog=logfolder+'/SEN_coreg_ifg_proc_stdout.log'

if Nburst == 0:
    graphxml=GRAPH+'/SEN_coreg_ifg_computation_subset_SINGLE_BURST.xml'
    print(graphxml)
    graph2run=PROJECT+'/graphs/SEN_coreg_ifg2run.xml'

    out_file = open(outlog, 'a')
    err_file=out_file

    print(bar_message)
    out_file.write(bar_message)
    message='## Coregistration and Interferogram computation started:\n'
    print(message)
    out_file.write(message)
    print(bar_message) 
    out_file.write(bar_message)
    print("# SINGLE BURST PROCESSING (no ESD operator) \n")
    k=0
    for dimfile in glob.iglob(slavesplittedfolder + '/*/*'+IW+'.dim'):
        print(dimfile)
        k=k+1
        head, tail = os.path.split(os.path.join(slavesplittedfolder, dimfile))
        message='['+str(k)+'] Processing slave file :'+tail+'\n'
        print(message)
        out_file.write(message)
        head , tailm = os.path.split(MASTER)
        outputname=tailm[17:25]+'_'+tail[0:8]+'_'+IW+'.dim'
        with open(graphxml, 'r') as file :
           filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('MASTER',MASTER)
        filedata = filedata.replace('SLAVE', dimfile)
        filedata = filedata.replace('DEMNAME',DEMNAME)
        filedata = filedata.replace('DEMFILE',DEMFILE)
        filedata = filedata.replace('COREGDTMNAME',COREGDTMNAME)
        filedata = filedata.replace('COREGDTMFILE',COREGDTMFILE)
        filedata = filedata.replace('DEM_RESAMPLING',DEMRESAMPLING)
        filedata = filedata.replace('OUTPUTCOREGFOLDER',outputcoregfolder)
        filedata = filedata.replace('OUTPUTIFGFOLDER', outputifgfolder)
        filedata = filedata.replace('OUTPUTFILE',outputname)
        filedata = filedata.replace('POLYGON',polygon)
        # Write the file out again
        with open(graph2run, 'w') as file:
           file.write(filedata)
        args = [ GPT, graph2run, '-c', CACHE, '-q', CPU]
        # Launch the processing
        process = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        timeStarted = time.time()
        stdout = process.communicate()[0]
        print('SNAP STDOUT:{}'.format(stdout))
        timeDelta = time.time() - timeStarted                     # Get execution time.
        print(('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.'))
        out_file.write('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.\n')
        if process.returncode != 0 :
            message='Error computing with coregistration and interferogram generation of splitted slave '+str(dimfile)
            err_file.write(message+'\n')
        else:
            message='Coregistration and Interferogram computation for data '+str(tail)+' successfully completed.\n'
            print(message)
            out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)
    out_file.close()
else:
    graphxml=GRAPH+'/SEN_coreg_ifg_computation_subset_MULTI_BURST.xml'
    print(graphxml)
    graph2run=PROJECT+'/graphs/SEN_coreg_ifg2run.xml'

    out_file = open(outlog, 'a')
    err_file=out_file

    print(bar_message)
    out_file.write(bar_message)
    message='## Coregistration and Interferogram computation started:\n'
    print(message)
    out_file.write(message)
    print(bar_message) 
    out_file.write(bar_message)
    print("# MULTI BURST PROCESSING (with ESD operator) \n")
    k=0
    for dimfile in glob.iglob(slavesplittedfolder + '/*/*'+IW+'.dim'):
        print(dimfile)
        k=k+1
        head, tail = os.path.split(os.path.join(slavesplittedfolder, dimfile))
        message='['+str(k)+'] Processing slave file :'+tail+'\n'
        print(message)
        out_file.write(message)
        head , tailm = os.path.split(MASTER)
        outputname=tailm[17:25]+'_'+tail[0:8]+'_'+IW+'.dim'
        with open(graphxml, 'r') as file :
           filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('MASTER',MASTER)
        filedata = filedata.replace('SLAVE', dimfile)
        filedata = filedata.replace('DEMNAME',DEMNAME)
        filedata = filedata.replace('DEMFILE',DEMFILE)
        filedata = filedata.replace('COREGDTMNAME',COREGDTMNAME)
        filedata = filedata.replace('COREGDTMFILE',COREGDTMFILE)
        filedata = filedata.replace('DEM_RESAMPLING',DEMRESAMPLING)
        filedata = filedata.replace('OUTPUTCOREGFOLDER',outputcoregfolder)
        filedata = filedata.replace('OUTPUTIFGFOLDER', outputifgfolder)
        filedata = filedata.replace('OUTPUTFILE',outputname)
        filedata = filedata.replace('POLYGON',polygon)
        # Write the file out again
        with open(graph2run, 'w') as file:
           file.write(filedata)
        args = [ GPT, graph2run, '-c', CACHE, '-q', CPU]
        # Launch the processing
        process = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        timeStarted = time.time()
        stdout = process.communicate()[0]
        print('SNAP STDOUT:{}'.format(stdout))
        timeDelta = time.time() - timeStarted                     # Get execution time.
        print(('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.'))
        out_file.write('['+str(k)+'] Finished process in '+str(timeDelta)+' seconds.\n')
        if process.returncode != 0 :
            message='Error computing with coregistration and interferogram generation of splitted slave '+str(dimfile)
            err_file.write(message+'\n')
        else:
            message='Coregistration and Interferogram computation for data '+str(tail)+' successfully completed.\n'
            print(message)
            out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)
    out_file.close()
