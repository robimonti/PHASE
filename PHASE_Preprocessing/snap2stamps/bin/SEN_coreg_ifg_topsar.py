### Python script to use SNAP as InSAR processor compatible with StaMPS PSI processing
# Author: Jose Manuel Delgado Blasco
# Modified by: Roberto Monti
# Automated Multi-Swath Dynamic Graph Generator (Memory-Safe Chunking & Band Preservation Edition)

import os
import sys
import glob
import subprocess
import time
import re
import shutil

inputfile = sys.argv[1]
bar_message='\n#####################################################################\n'

print("############################# STEP 3 ################################\n")
print("## INPUT PARAMETERS ##")
try:
    in_file = open(inputfile, 'r')
    for line in in_file.readlines():
        if "PROJECTFOLDER" in line: PROJECT = line.split('=')[1].strip()
        if "DEMNAME" in line: DEMNAME = line.split('=')[1].strip()
        if "DEMFILE" in line: DEMFILE = line.split('=')[1].strip()
        if "COREGDTMNAME" in line: COREGDTMNAME = line.split('=')[1].strip()
        if "COREGDTMFILE" in line: COREGDTMFILE = line.split('=')[1].strip()
        if "DEMRESAMPLING" in line: DEMRESAMPLING = line.split('=')[1].strip()
        if "SWATHS" in line:
            SWATHS_STR = line.split('=')[1].strip()
            swaths = [s.strip() for s in SWATHS_STR.split(',')]
        if "MASTER" in line: MASTER = line.split('=')[1].strip()
        if "GRAPHSFOLDER" in line: GRAPH = line.split('=')[1].strip()
        if "GPTBIN_PATH" in line: GPT = line.split('=')[1].strip()
        if "LONMIN" in line: LONMIN = line.split('=')[1].strip()
        if "LATMIN" in line: LATMIN = line.split('=')[1].strip()
        if "LONMAX" in line: LONMAX = line.split('=')[1].strip()
        if "LATMAX" in line: LATMAX = line.split('=')[1].strip()
        if "CACHE" in line: CACHE = line.split('=')[1].strip()
        if "CPU" in line: CPU = line.split('=')[1].strip()
finally:
    in_file.close()

# 1. Geometry and Path Extraction
polygon = f"POLYGON (({LONMIN} {LATMIN},{LONMAX} {LATMIN},{LONMAX} {LATMAX},{LONMIN} {LATMAX},{LONMIN} {LATMIN}))"

mastersplittedfolder = os.path.dirname(MASTER)
master_filename = os.path.basename(MASTER)

# Split the filename at '_split_' to isolate the 67-character prefix
master_prefix = master_filename.split('_split_')[0]
master_date = master_prefix[17:25]

slavesplittedfolder = os.path.join(PROJECT, 'split')
outputcoregfolder = os.path.join(PROJECT, 'coreg')
outputifgfolder = os.path.join(PROJECT, 'ifg')
logfolder = os.path.join(PROJECT, 'logs')
graphsfolder = os.path.join(PROJECT, 'graphs')

for f in [outputcoregfolder, outputifgfolder, logfolder, graphsfolder]:
    if not os.path.exists(f):
        os.makedirs(f)

outlog = os.path.join(logfolder, 'SEN_coreg_ifg_proc_stdout.log')
out_file = open(outlog, 'a')
err_file = out_file

def count_bursts(dim_path):
    """Aggressively parses the .dim XML to find the true burst count, reading the literal burst tags."""
    try:
        with open(dim_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # The Ground Truth: Count the literal <MDElem name="burst"> tags
        bursts = len(re.findall(r'<MDElem name="burst">', content))
        if bursts > 0:
            return bursts
        
        # Fallbacks for older or alternative SNAP formats
        match = re.search(r'name="numBursts"[^>]*>\s*(\d+)\s*<', content, re.IGNORECASE)
        if match: return int(match.group(1))
        
        match = re.search(r'<num_?bursts[^>]*>\s*(\d+)\s*</num_?bursts>', content, re.IGNORECASE)
        if match: return int(match.group(1))
        
        return 1
    except Exception as e:
        print(f"Warning: Burst count failed for {dim_path}. Defaulting to 1. Error: {e}")
        return 1

print(bar_message)
out_file.write(bar_message)
message='## Dynamic Coregistration and Interferogram computation started:\n'
print(message)
out_file.write(message)
print(bar_message)
out_file.write(bar_message)

k = 0
for slave_date in sorted(os.listdir(slavesplittedfolder)):
    slave_date_path = os.path.join(slavesplittedfolder, slave_date)
    if not os.path.isdir(slave_date_path):
        continue

    valid_swaths = []
    for swath in swaths:
        master_dim = os.path.join(mastersplittedfolder, f"{master_prefix}_split_{swath}_Orb.dim")
        slave_dim = os.path.join(slave_date_path, f"{slave_date}_{swath}.dim")
        if os.path.exists(master_dim) and os.path.exists(slave_dim):
            valid_swaths.append(swath)

    if not valid_swaths:
        continue 

    k += 1
    print(f"\n[{k}] Processing Master {master_date} -> Slave {slave_date}")
    print(f"    Valid intersecting swaths detected: {valid_swaths}")
    out_file.write(f"\n[{k}] Processing Master {master_date} -> Slave {slave_date}\n")

    outputname = f"{master_date}_{slave_date}.dim"
    coreg_out = os.path.join(outputcoregfolder, outputname)
    ifg_out = os.path.join(outputifgfolder, outputname)

    temp_coregs = []
    temp_ifgs = []

    # --- PHASE 1: PROCESS EACH SWATH SEQUENTIALLY ---
    for swath in valid_swaths:
        master_dim = os.path.join(mastersplittedfolder, f"{master_prefix}_split_{swath}_Orb.dim")
        slave_dim = os.path.join(slave_date_path, f"{slave_date}_{swath}.dim")
        
        n_bursts = count_bursts(master_dim)
        
        xml = ['<graph id="Graph">\n  <version>1.0</version>']
        xml.append(f'''
  <node id="Read-Master-{swath}">
    <operator>Read</operator>
    <parameters><file>{master_dim}</file></parameters>
  </node>
  <node id="Read-Slave-{swath}">
    <operator>Read</operator>
    <parameters><file>{slave_dim}</file></parameters>
  </node>
  <node id="Back-Geocoding-{swath}">
    <operator>Back-Geocoding</operator>
    <sources>
      <sourceProduct refid="Read-Master-{swath}"/>
      <sourceProduct.1 refid="Read-Slave-{swath}"/>
    </sources>
    <parameters>
      <demName>{DEMNAME}</demName>
      <demResamplingMethod>{DEMRESAMPLING}</demResamplingMethod>
      <externalDEMFile>{DEMFILE}</externalDEMFile>
      <externalDEMNoDataValue>0.0</externalDEMNoDataValue>
      <resamplingType>BILINEAR_INTERPOLATION</resamplingType>
    </parameters>
  </node>''')

        last_coreg_node = f"Back-Geocoding-{swath}"

        if n_bursts > 1:
            print(f"    -> {swath}: {n_bursts} bursts detected. Injecting ESD Node.")
            xml.append(f'''
  <node id="ESD-{swath}">
    <operator>Enhanced-Spectral-Diversity</operator>
    <sources>
      <sourceProduct refid="{last_coreg_node}"/>
    </sources>
    <parameters>
      <fineWinWidthStr>512</fineWinWidthStr>
      <fineWinHeightStr>512</fineWinHeightStr>
    </parameters>
  </node>''')
            last_coreg_node = f"ESD-{swath}"
        else:
            print(f"    -> {swath}: 1 burst detected. Bypassing ESD Node.")

        # REMOVED the restrictive <parameters> blocks from Deburst to allow elevation to pass
        xml.append(f'''
  <node id="Interferogram-{swath}">
    <operator>Interferogram</operator>
    <sources>
      <sourceProduct refid="{last_coreg_node}"/>
    </sources>
    <parameters>
      <subtractFlatEarthPhase>true</subtractFlatEarthPhase>
      <srpPolynomialDegree>5</srpPolynomialDegree>
      <srpNumberPoints>501</srpNumberPoints>
      <orbitDegree>3</orbitDegree>
      <includeCoherence>true</includeCoherence>
      <outputElevation>true</outputElevation>
      <outputLatLon>true</outputLatLon>
      <demName>{DEMNAME}</demName>
      <externalDEMFile>{DEMFILE}</externalDEMFile>
      <externalDEMNoDataValue>0.0</externalDEMNoDataValue>
    </parameters>
  </node>
  <node id="Deburst-Coreg-{swath}">
    <operator>TOPSAR-Deburst</operator>
    <sources>
      <sourceProduct refid="{last_coreg_node}"/>
    </sources>
  </node>
  <node id="Deburst-IFG-{swath}">
    <operator>TOPSAR-Deburst</operator>
    <sources>
      <sourceProduct refid="Interferogram-{swath}"/>
    </sources>
  </node>''')

        temp_c_name = os.path.join(outputcoregfolder, f"temp_{master_date}_{slave_date}_coreg_{swath}.dim")
        temp_i_name = os.path.join(outputifgfolder, f"temp_{master_date}_{slave_date}_ifg_{swath}.dim")
        temp_coregs.append(temp_c_name)
        temp_ifgs.append(temp_i_name)

        xml.append(f'''
  <node id="Write-Temp-Coreg-{swath}">
    <operator>Write</operator>
    <sources><sourceProduct refid="Deburst-Coreg-{swath}"/></sources>
    <parameters><file>{temp_c_name}</file><formatName>BEAM-DIMAP</formatName></parameters>
  </node>
  <node id="Write-Temp-IFG-{swath}">
    <operator>Write</operator>
    <sources><sourceProduct refid="Deburst-IFG-{swath}"/></sources>
    <parameters><file>{temp_i_name}</file><formatName>BEAM-DIMAP</formatName></parameters>
  </node>
</graph>''')

        graph_swath = os.path.join(graphsfolder, f"SEN_temp_{swath}.xml")
        with open(graph_swath, 'w') as f: f.write("\n".join(xml))
        
        print(f"    -> Processing {swath} ...")
        args = [GPT, graph_swath, '-c', CACHE, '-q', CPU]
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout = process.communicate()[0].decode('utf-8', errors='ignore')
        if process.returncode != 0:
            print(f"ERROR processing {swath}:\n{stdout}")
            sys.exit(1)

    # --- PHASE 2: MERGE, TOPOPHASE, SUBSET, FINAL WRITE ---
    xml = ['<graph id="Graph">\n  <version>1.0</version>']
    
    coreg_merge_sources = ""
    ifg_merge_sources = ""
    
    for i, (tc, ti) in enumerate(zip(temp_coregs, temp_ifgs)):
        xml.append(f'''
  <node id="Read-Coreg-{i}">
    <operator>Read</operator>
    <parameters><file>{tc}</file></parameters>
  </node>
  <node id="Read-IFG-{i}">
    <operator>Read</operator>
    <parameters><file>{ti}</file></parameters>
  </node>''')
        src_prefix = "sourceProduct" if i == 0 else f"sourceProduct.{i}"
        coreg_merge_sources += f'      <{src_prefix} refid="Read-Coreg-{i}"/>\n'
        ifg_merge_sources += f'      <{src_prefix} refid="Read-IFG-{i}"/>\n'

    if len(valid_swaths) > 1:
        print("    -> Generating TOPSAR-Merge nodes for combined footprint.")
        # REMOVED parameters block to prevent Merge from deleting topo bands
        xml.append(f'''
  <node id="Merge-Coreg">
    <operator>TOPSAR-Merge</operator>
    <sources>\n{coreg_merge_sources}    </sources>
  </node>
  <node id="Merge-IFG">
    <operator>TOPSAR-Merge</operator>
    <sources>\n{ifg_merge_sources}    </sources>
  </node>''')
        final_coreg_node = "Merge-Coreg"
        final_ifg_node = "Merge-IFG"
    else:
        final_coreg_node = "Read-Coreg-0"
        final_ifg_node = "Read-IFG-0"

    # BELT AND SUSPENDERS: Forced TopoPhaseRemoval to explicitly output elevation/lat/lon
    xml.append(f'''
  <node id="TopoPhaseRemoval">
    <operator>TopoPhaseRemoval</operator>
    <sources>
      <sourceProduct refid="{final_ifg_node}"/>
    </sources>
    <parameters>
      <orbitDegree>3</orbitDegree>
      <demName>{DEMNAME}</demName>
      <externalDEMFile>{DEMFILE}</externalDEMFile>
      <externalDEMNoDataValue>0.0</externalDEMNoDataValue>
      <outputTopoPhaseBand>false</outputTopoPhaseBand>
      <outputElevationBand>true</outputElevationBand>
      <outputLatLonBands>true</outputLatLonBands>
    </parameters>
  </node>
  <node id="Subset-Coreg">
    <operator>Subset</operator>
    <sources>
      <sourceProduct refid="{final_coreg_node}"/>
    </sources>
    <parameters>
      <geoRegion>{polygon}</geoRegion>
      <copyMetadata>true</copyMetadata>
    </parameters>
  </node>
  <node id="Subset-IFG">
    <operator>Subset</operator>
    <sources>
      <sourceProduct refid="TopoPhaseRemoval"/>
    </sources>
    <parameters>
      <geoRegion>{polygon}</geoRegion>
      <copyMetadata>true</copyMetadata>
    </parameters>
  </node>
  <node id="Write-Coreg">
    <operator>Write</operator>
    <sources>
      <sourceProduct refid="Subset-Coreg"/>
    </sources>
    <parameters>
      <file>{coreg_out}</file>
      <formatName>BEAM-DIMAP</formatName>
    </parameters>
  </node>
  <node id="Write-IFG">
    <operator>Write</operator>
    <sources>
      <sourceProduct refid="Subset-IFG"/>
    </sources>
    <parameters>
      <file>{ifg_out}</file>
      <formatName>BEAM-DIMAP</formatName>
    </parameters>
  </node>
</graph>''')

    print("    -> Finalizing Merge, TopoPhaseRemoval, and Subsetting...")
    graph_final = os.path.join(graphsfolder, f"SEN_final_{slave_date}.xml")
    with open(graph_final, 'w') as f: f.write("\n".join(xml))
    
    args = [GPT, graph_final, '-c', CACHE, '-q', CPU]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    timeStarted = time.time()
    stdout = process.communicate()[0].decode('utf-8', errors='ignore')
    timeDelta = time.time() - timeStarted

    if process.returncode != 0:
        msg = f'    -> ERROR computing final graph for {slave_date}\n{stdout}\n'
        err_file.write(msg)
        print(msg)
    else:
        msg = f'    -> Coregistration completed in {timeDelta:.2f} seconds.\n'
        print(msg)
        out_file.write(msg)
        
        # Clean up temporary disk files
        for tc, ti in zip(temp_coregs, temp_ifgs):
            os.remove(tc)
            shutil.rmtree(tc.replace('.dim', '.data'))
            os.remove(ti)
            shutil.rmtree(ti.replace('.dim', '.data'))

out_file.close()