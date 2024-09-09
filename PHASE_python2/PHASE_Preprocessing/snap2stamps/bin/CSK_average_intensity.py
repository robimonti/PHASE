import os
from pathlib import Path
import sys
import subprocess
import shlex
import time

inputfile = sys.argv[1]

bar_message = '\n#####################################################################\n'

# Getting configuration variables from inputfile
print("############################# STEP 5 ################################\n")
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
        if "CACHE" in line:
            CACHE = line.split('=')[1].strip()
        if "CPU" in line:
            CPU = line.split('=')[1].strip()
finally:
    in_file.close()

#############################################################################
############@##### Average intensity image and export (tif) #################
#############################################################################
intensityfolder = PROJECT + '/intensity'
subsetfolder = PROJECT + '/subset'
masterfolder = PROJECT + '/master'
logfolder = PROJECT + '/logs'
graphfolder = PROJECT + '/graphs'
graphs_dir = os.path.join(PROJECT, 'snap2stamps', 'graphs')
if not os.path.exists(intensityfolder):
    os.makedirs(intensityfolder)

outlog = logfolder + '/CSK_avg_intensity_proc_stdout.log'

graphxml = GRAPH + '/CSK_slave_intensity.xml'
print(graphxml)
graph2run = PROJECT + '/graphs/CSK_slave_intensity2run.xml'

# Find all slave folders
slave_folders = [os.path.join(subsetfolder, folder) for folder in os.listdir(subsetfolder) if os.path.isdir(os.path.join(subsetfolder, folder))]

# Find all .dim files in each slave folder
slave_files = []
for folder in slave_folders:
    dim_slaves = [os.path.join(folder, filename) for filename in os.listdir(folder) if filename.endswith(".dim")]
    slave_files.extend(dim_slaves)

# Find master image
for filename in os.listdir(masterfolder):
    if filename.endswith("_sub.dim"):
        master_product = os.path.join(masterfolder, filename)
        break

# Ensure the master product is included as the first product in the stack
slave_files.insert(0, master_product)

# Print the list of slave files (including the master as the first item)
for file in slave_files:
    print(file)
    
# Count the number of files
Nstack = len(slave_files)


# Create an XML string
xml_content = '''
<graph id="Graph">
  <version>1.0</version>
'''

# Create Nstack Read nodes and connect them to the CreateStack node
for index, file in enumerate(slave_files):
    xml_content += '''
  <node id="Read_{}">
    <operator>Read</operator>
    <sources/>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>{}</file>
    </parameters>
  </node>
    '''.format(index, file)

xml_content += '''
  <node id="CreateStack">
    <operator>CreateStack</operator>
    <sources>
      <sourceProduct refid="Read_0"/>
'''

# Connect the Read nodes to the CreateStack node with .n for n > 0
for index in range(1, Nstack):
    xml_content += '''
      <sourceProduct.{} refid="Read_{}"/>
    '''.format(index, index)

xml_content += '''
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <masterBands/>
      <sourceBands/>
      <resamplingType>NONE</resamplingType>
      <extent>Master</extent>
      <initialOffsetMethod>Orbit</initialOffsetMethod>
    </parameters>
  </node>
'''

# Add the Stack Averaging node
xml_content += '''
  <node id="Stack-Averaging">
    <operator>Stack-Averaging</operator>
    <sources>
      <sourceProduct refid="CreateStack"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <statistic>Mean Average</statistic>
    </parameters>
  </node>
'''

# Add the Band Maths node
xml_content += '''
  <node id="BandMaths">
    <operator>BandMaths</operator>
    <sources>
      <sourceProduct refid="Stack-Averaging"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <targetBands>
        <targetBand>
          <name>avgIntensity</name>
          <type>float32</type>
          <expression>i==0.0?0.0:i*i+q*q</expression>
          <description/>
          <unit/>
          <noDataValue>0.0</noDataValue>
        </targetBand>
      </targetBands>
      <variables/>
    </parameters>
  </node>
'''

# Add the Terrain Correction node
xml_content += '''
  <node id="Terrain-Correction">
    <operator>Terrain-Correction</operator>
    <sources>
      <sourceProduct refid="BandMaths"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <sourceBands>avgIntensity</sourceBands>
      <demName>DEMNAME</demName>
      <externalDEMFile>DEMFILE</externalDEMFile>
      <externalDEMNoDataValue>0.0</externalDEMNoDataValue>
      <externalDEMApplyEGM>false</externalDEMApplyEGM>
      <demResamplingMethod>BILINEAR_INTERPOLATION</demResamplingMethod>
      <imgResamplingMethod>BILINEAR_INTERPOLATION</imgResamplingMethod>
      <pixelSpacingInMeter>2.2058094909171118</pixelSpacingInMeter>
      <pixelSpacingInDegree>1.9815123795467424E-5</pixelSpacingInDegree>
      <mapProjection>GEOGCS["WGS84(DD)", 
			DATUM["WGS84", 
				SPHEROID["WGS84", 6378137.0, 298.257223563]], 
			PRIMEM["Greenwich", 0.0], 
			UNIT["degree", 0.017453292519943295], 
			AXIS["Geodetic longitude", EAST], 
			AXIS["Geodetic latitude", NORTH]]</mapProjection>
      <alignToStandardGrid>false</alignToStandardGrid>
      <standardGridOriginX>0.0</standardGridOriginX>
      <standardGridOriginY>0.0</standardGridOriginY>
      <nodataValueAtSea>true</nodataValueAtSea>
      <saveDEM>false</saveDEM>
      <saveLatLon>false</saveLatLon>
      <saveIncidenceAngleFromEllipsoid>false</saveIncidenceAngleFromEllipsoid>
      <saveLocalIncidenceAngle>false</saveLocalIncidenceAngle>
      <saveProjectedLocalIncidenceAngle>false</saveProjectedLocalIncidenceAngle>
      <saveSelectedSourceBand>true</saveSelectedSourceBand>
      <saveLayoverShadowMask>false</saveLayoverShadowMask>
      <outputComplex>false</outputComplex>
      <applyRadiometricNormalization>false</applyRadiometricNormalization>
      <saveSigmaNought>false</saveSigmaNought>
      <saveGammaNought>false</saveGammaNought>
      <saveBetaNought>false</saveBetaNought>
      <incidenceAngleForSigma0>Use projected local incidence angle from DEM</incidenceAngleForSigma0>
      <incidenceAngleForGamma0>Use projected local incidence angle from DEM</incidenceAngleForGamma0>
      <auxFile>Latest Auxiliary File</auxFile>
      <externalAuxFile/>
    </parameters>
  </node>
'''

# Add the Write node
xml_content += '''
  <node id="Write">
    <operator>Write</operator>
    <sources>
      <sourceProduct refid="BandMaths"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>OUTPUTINTFOLDER/OUTPUTFILE_INT</file>
      <formatName>BEAM-DIMAP</formatName>
    </parameters>
  </node>
'''

# Add the Write GEOTIFF node
xml_content += '''
  <node id="Write1">
    <operator>Write</operator>
    <sources>
      <sourceProduct refid="Terrain-Correction"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <file>OUTPUTINTFOLDER/OUTPUTFILE_INT_TC</file>
      <formatName>GeoTIFF</formatName>
    </parameters>
  </node>
'''

# Complete the XML content
xml_content += '</graph>'

xml_file_path = os.path.join(graphs_dir, 'CSK_slave_intensity.xml')
with open(xml_file_path, 'w') as xml_file:
    xml_file.write(xml_content)

print('Generated XML saved as CSK_slave_intensity.xml')


print(bar_message)
with open(outlog, 'a') as out_file:
    err_file = out_file

    out_file.write(bar_message)
    message = '## Average Intensity scene \n'
    print(message)
    out_file.write(message)
    print(bar_message)
    out_file.write(bar_message)

    
    outputname_int = 'AverageIntensity'
    outputname_int_tc = 'AverageIntensity_TC.tif'
    with open(graphxml, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('{{SLAVES}}', ','.join(slave_files))
    filedata = filedata.replace('DEMNAME', DEMNAME)
    filedata = filedata.replace('DEMFILE', DEMFILE)
    filedata = filedata.replace('OUTPUTFILE_INT_TC', outputname_int_tc)
    filedata = filedata.replace('OUTPUTFILE_INT', outputname_int)
    filedata = filedata.replace('OUTPUTINTFOLDER', intensityfolder)

    # Write the file out again
    with open(graph2run, 'w') as file:
        file.write(filedata)

    args = [GPT, graph2run, '-c', CACHE, '-q', CPU]

    # Launch the processing
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    timeStarted = time.time()
    stdout = process.communicate()[0]
    print('SNAP STDOUT:{}'.format(stdout))
    timeDelta = time.time() - timeStarted  # Get execution time.
    print('Finished process in ' + str(timeDelta) + ' seconds.')
    out_file.write('Finished process in ' + str(timeDelta) + ' seconds.\n')

    if process.returncode != 0:
        message = 'Error average intensity '
        err_file.write(message)
    else:
        message = 'Average intensity image successfully completed.\n'
        print(message)
        out_file.write(message)
        print(bar_message)
        out_file.write(bar_message)

