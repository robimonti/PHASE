### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Preserves Date Subfolder Architecture, upgraded with Regex Date Extraction

import os
import shutil
import sys
import re

inputfile = sys.argv[1]
bar_message='\n#####################################################################\n'

print("############################# STEP 1 ################################\n")
print("## INPUT PARAMETERS ##")
try:
    with open(inputfile, 'r') as in_file:
        for line in in_file.readlines():
            if "PROJECTFOLDER" in line:
                PROJECT = line.split('=')[1].strip()
                print("\n- PROJECT FOLDER: ")
                print(PROJECT)
except Exception as e:
    print(f"Error reading config: {e}")
    sys.exit(1)

logfolder = os.path.join(PROJECT, 'logs')
if not os.path.exists(logfolder):
    os.makedirs(logfolder)

outlog = os.path.join(logfolder, 'split_proc_stdout.log')

with open(outlog, 'a') as out_file:
    print(bar_message)
    out_file.write(bar_message)
    message = '## Slaves sorting into folders\n'
    print(message)
    out_file.write(message)

    directory = os.path.join(PROJECT, 'slaves')
    for filename in os.listdir(directory):
        if filename.endswith(".h5"):
            # Safely extract Date using Regex instead of [37:45]
            match = re.search(r'_(\d{8})\d{6}_', filename)
            date_str = match.group(1) if match else filename[37:45]

            print(f"Processing: {filename} -> Date: {date_str}")

            subdirectory = os.path.join(directory, date_str)
            if not os.path.exists(subdirectory):
                os.makedirs(subdirectory)

            source = os.path.join(directory, filename)
            destination = os.path.join(subdirectory, filename)
            print(f'Moving {source} to {destination}')
            shutil.move(source, destination)
