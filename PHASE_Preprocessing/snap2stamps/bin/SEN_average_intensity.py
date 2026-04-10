### Python script to compute Full-Stack Average Intensity for Sentinel-1
# Upgraded to bypass SNAP 9 BandMerge conflicts and robust Band Extraction

import os
import sys
import subprocess
import time
import glob
import re

inputfile = sys.argv[1]
bar_message = '\n#####################################################################\n'

print("############################# STEP 5 ################################\n")
print("## INPUT PARAMETERS ##")
try:
    with open(inputfile, 'r') as in_file:
        for line in in_file.readlines():
            if "PROJECTFOLDER" in line: PROJECT = line.split('=')[1].strip()
            if "DEMNAME" in line: DEMNAME = line.split('=')[1].strip()
            if "DEMFILE" in line: DEMFILE = line.split('=')[1].strip()
            if "DEMRESAMPLING" in line: DEMRESAMPLING = line.split('=')[1].strip()
            if "LONMIN" in line: LONMIN = line.split('=')[1].strip()
            if "LATMIN" in line: LATMIN = line.split('=')[1].strip()
            if "LONMAX" in line: LONMAX = line.split('=')[1].strip()
            if "LATMAX" in line: LATMAX = line.split('=')[1].strip()
            if "GPTBIN_PATH" in line: GPT = line.split('=')[1].strip()
            if "CACHE" in line: CACHE = line.split('=')[1].strip()
            if "CPU" in line: CPU = line.split('=')[1].strip()
except Exception as e:
    print(f"Error reading config: {e}")
    sys.exit(1)

polygon = f"POLYGON (({LONMIN} {LATMIN},{LONMAX} {LATMIN},{LONMAX} {LATMAX},{LONMIN} {LATMAX},{LONMIN} {LATMIN}))"

intensityfolder = os.path.join(PROJECT, 'intensity')
coregfolder = os.path.join(PROJECT, 'coreg')
logfolder = os.path.join(PROJECT, 'logs')
graphs_dir = os.path.join(PROJECT, 'graphs')

for f in [intensityfolder, logfolder, graphs_dir]:
    if not os.path.exists(f):
        os.makedirs(f)

outlog = os.path.join(logfolder, 'SEN_avg_intensity_proc_stdout.log')
graphxml = os.path.join(graphs_dir, 'SEN_stack_intensity.xml')

# 1. Grab all Coregistered files
coreg_files = glob.glob(os.path.join(coregfolder, "*.dim"))
if not coreg_files:
    print("ERROR: No coregistered .dim files found. Run Step 3 first.")
    sys.exit(1)

print(f"- Found {len(coreg_files)} coregistered pairs for Full-Stack Average Intensity.")

# 2. BULLETPROOF DYNAMIC BAND EXTRACTION
def get_iq_bands(dim_file):
    with open(dim_file, 'r', encoding='utf-8') as f:
        content = f.read()
    bands = re.findall(r'<BAND_NAME>(.*?)</BAND_NAME>', content)
    # Simple startswith bypasses any missing 'VV' or 'IW' tags
    i_b = [b for b in bands if b.startswith('i_')]
    q_b = [b for b in bands if b.startswith('q_')]
    return i_b, q_b

i_bands_all = []
q_bands_all = []
for f in coreg_files:
    i_b, q_b = get_iq_bands(f)
    i_bands_all.append(i_b)
    q_bands_all.append(q_b)

# Find Master bands mathematically
master_i_set = set(i_bands_all[0])
master_q_set = set(q_bands_all[0])

if len(coreg_files) > 1:
    for i in range(1, len(i_bands_all)):
        master_i_set = master_i_set.intersection(i_bands_all[i])
        master_q_set = master_q_set.intersection(q_bands_all[i])
        
if len(master_i_set) == 1:
    master_i_name = list(master_i_set)[0]
    master_q_name = list(master_q_set)[0]
else:
    # Fallback to prevent IndexError
    master_i_name = [b for b in i_bands_all[0] if 'slv' not in b.lower() and 'late' not in b.lower()][0]
    master_q_name = [b for b in q_bands_all[0] if 'slv' not in b.lower() and 'late' not in b.lower()][0]

print(f"- Dynamically detected Master Bands: {master_i_name}, {master_q_name}")

# 3. Build the Pre-Calculated XML Graph
xml_content = '<graph id="Graph">\n  <version>1.0</version>\n'

for index, file in enumerate(coreg_files):
    xml_content += f'''
  <node id="Read_{index}">
    <operator>Read</operator>
    <parameters><file>{file}</file></parameters>
  </node>'''

xml_content += f'''
  <node id="Math_Master">
    <operator>BandMaths</operator>
    <sources><sourceProduct refid="Read_0"/></sources>
    <parameters>
      <targetBands>
        <targetBand>
          <name>Int_Master</name>
          <type>float32</type>
          <expression>{master_i_name} * {master_i_name} + {master_q_name} * {master_q_name}</expression>
        </targetBand>
      </targetBands>
    </parameters>
  </node>'''

merge_sources = ['      <sourceProduct refid="Math_Master"/>']
avg_expression_parts = ['Int_Master']

for index, file in enumerate(coreg_files):
    slave_i_name = [b for b in i_bands_all[index] if b != master_i_name][0]
    slave_q_name = [b for b in q_bands_all[index] if b != master_q_name][0]
    
    xml_content += f'''
  <node id="Math_Slave_{index}">
    <operator>BandMaths</operator>
    <sources><sourceProduct refid="Read_{index}"/></sources>
    <parameters>
      <targetBands>
        <targetBand>
          <name>Int_Slave_{index}</name>
          <type>float32</type>
          <expression>{slave_i_name} * {slave_i_name} + {slave_q_name} * {slave_q_name}</expression>
        </targetBand>
      </targetBands>
    </parameters>
  </node>'''
    
    merge_sources.append(f'      <sourceProduct.{index+1} refid="Math_Slave_{index}"/>')
    avg_expression_parts.append(f'Int_Slave_{index}')

merge_sources_str = "\n".join(merge_sources)
xml_content += f'''
  <node id="BandMerge">
    <operator>BandMerge</operator>
    <sources>
{merge_sources_str}
    </sources>
    <parameters><geographicError>1.0E-5</geographicError></parameters>
  </node>'''

total_images = len(coreg_files) + 1
expr_str = " + ".join(avg_expression_parts)
final_expr = f'({expr_str}) / {total_images}'

xml_content += f'''
  <node id="Math_Avg">
    <operator>BandMaths</operator>
    <sources><sourceProduct refid="BandMerge"/></sources>
    <parameters>
      <targetBands>
        <targetBand>
          <name>avgIntensity</name>
          <type>float32</type>
          <expression>{final_expr}</expression>
        </targetBand>
      </targetBands>
    </parameters>
  </node>'''

outputname_int = os.path.join(intensityfolder, 'AverageIntensity.dim')
outputname_int_tc = os.path.join(intensityfolder, 'AverageIntensity_TC.tif')

# S1 Needs the Subset node!
xml_content += f'''
  <node id="Subset">
    <operator>Subset</operator>
    <sources><sourceProduct refid="Math_Avg"/></sources>
    <parameters>
      <geoRegion>{polygon}</geoRegion>
      <copyMetadata>true</copyMetadata>
    </parameters>
  </node>
  <node id="Terrain-Correction">
    <operator>Terrain-Correction</operator>
    <sources><sourceProduct refid="Subset"/></sources>
    <parameters>
      <sourceBands>avgIntensity</sourceBands>
      <demName>{DEMNAME}</demName>
      <externalDEMFile>{DEMFILE}</externalDEMFile>
      <externalDEMNoDataValue>0.0</externalDEMNoDataValue>
      <demResamplingMethod>{DEMRESAMPLING}</demResamplingMethod>
      <imgResamplingMethod>{DEMRESAMPLING}</imgResamplingMethod>
      <pixelSpacingInMeter>13.95755</pixelSpacingInMeter>
      <saveSelectedSourceBand>true</saveSelectedSourceBand>
    </parameters>
  </node>
  <node id="Write-DIM">
    <operator>Write</operator>
    <sources><sourceProduct refid="Subset"/></sources>
    <parameters><file>{outputname_int}</file><formatName>BEAM-DIMAP</formatName></parameters>
  </node>
  <node id="Write-TIF">
    <operator>Write</operator>
    <sources><sourceProduct refid="Terrain-Correction"/></sources>
    <parameters><file>{outputname_int_tc}</file><formatName>GeoTIFF</formatName></parameters>
  </node>
</graph>'''

with open(graphxml, 'w') as xml_file:
    xml_file.write(xml_content)

print(bar_message)
print("## Computing Full-Stack Average Intensity scene")
print(f"-> Merging Master + {len(coreg_files)} Slaves...")

with open(outlog, 'a') as out_file:
    out_file.write("## Average Intensity Processing\n")
    args = [GPT, graphxml, '-c', CACHE, '-q', CPU]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    timeStarted = time.time()
    stdout = process.communicate()[0].decode('utf-8', errors='ignore')
    timeDelta = time.time() - timeStarted
    
    if process.returncode != 0:
        print("-> ERROR: Average intensity computation failed.")
        print(f'SNAP STDOUT: {stdout}')
        out_file.write(f"ERROR. Execution time: {timeDelta:.2f}s\n")
    else:
        print(f"-> Success: Average intensity image completed in {timeDelta:.2f} seconds.")
        out_file.write(f"Success. Execution time: {timeDelta:.2f}s\n")

print(bar_message)