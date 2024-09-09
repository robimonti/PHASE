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
        if "IW1" in line:
            IWs = line.split('=')[1].strip()
            print("\n- SWATH: ")
            print IWs
        if "POLARISATION" in line:
            POLARISATION = line.split('=')[1].strip()
            print("\n- POLARIZATION: ")
            print POLARISATION
        if "DEMNAME" in line:
            DEMNAME = line.split('=')[1].strip()
            print("\n- DEM NAME: ")
            print(DEMNAME)
        if "DEMFILE" in line:
            DEMFILE = line.split('=')[1].strip()
            print("\n- DEM FILE: ")
            print(DEMFILE)
        if "DEMRESAMPLING" in line:
            DEMRESAMPLING = line.split('=')[1].strip()
            print("\n- DEM RESAMPLING METHOD: ")
            print(DEMRESAMPLING)
        if "GRAPHSFOLDER" in line:
            GRAPH = line.split('=')[1].strip()
            print("\n- GRAPH FOLDER: ")
            print(GRAPH)
        if "LONMIN" in line:
			LONMIN = line.split('=')[1].strip()
        if "LATMIN" in line:
            LATMIN = line.split('=')[1].strip()
        if "LONMAX" in line:
            LONMAX = line.split('=')[1].strip()
        if "LATMAX" in line:
            LATMAX = line.split('=')[1].strip()
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

polygon='POLYGON (('+LONMIN+' '+LATMIN+','+LONMAX+' '+LATMIN+','+LONMAX+' '+LATMAX+','+LONMIN+' '+LATMAX+','+LONMIN+' '+LATMIN+'))'

#############################################################################
############@##### Average intensity image and export (tif) #################
#############################################################################
intensityfolder = PROJECT + '/intensity'
coregfolder = PROJECT + '/coreg'
masterfolder = PROJECT + '/master'
logfolder = PROJECT + '/logs'
graphfolder = PROJECT + '/graphs'
graphs_dir = os.path.join(PROJECT, 'snap2stamps', 'graphs')
if not os.path.exists(intensityfolder):
    os.makedirs(intensityfolder)

outlog = logfolder + '/SEN_avg_intensity_proc_stdout.log'

graphxml = GRAPH + '/SEN_slave_intensity.xml'
print(graphxml)
graph2run = PROJECT + '/graphs/SEN_slave_intensity2run.xml'

# Check if the XML file exists and delete it if it does
if os.path.exists(graphxml):
    os.remove(graphxml)
if os.path.exists(graph2run):
    os.remove(graph2run)

# Find all slave files
slave_files = [os.path.join(coregfolder, filename) for filename in os.listdir(coregfolder) if filename.endswith(".dim")]

# Find master image
for filename in os.listdir(masterfolder):
    if filename.endswith("_split_Orb.dim"):
        master_product = os.path.join(masterfolder, filename)
        break

# Extract the date part from the master image filename
#master_filename = os.path.basename(master_product)
#master_date_part = master_filename.split("_")[7]

# Convert the date format for the master image (e.g., "20190721" to "21Jul2019")
#master_date_formatted = time.strptime(master_date_part, "%Y%m%d")
#master_date_formatted = time.strftime("%d%b%Y", master_date_formatted)

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

# Create Nstack Read nodes and BandsExtractorOp nodes
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
    
# Add the BandsExtractorOp node for the master image
xml_content += '''
  <node id="BandsExtractor_0">
    <operator>BandsExtractorOp</operator>
    <sources>
      <sourceProduct refid="Read_0"/>
    </sources>
    <parameters>
      <sourceBandNames>i_{}_{},q_{}_{}</sourceBandNames>
    </parameters>
  </node>
'''.format(IWs, POLARISATION, IWs, POLARISATION)

# Create BandsExtractorOp nodes for slave products
for index, file in enumerate(slave_files[1:], start=1):
    # Extract the date part from the filename in the file path
    filename = os.path.basename(file)  # Get just the filename
    date_part = filename.split("_")[1]  # Assuming the date is the second part between underscores
    # Convert the date format (e.g., "20190709" to "09Jul2019")
    date_formatted = time.strptime(date_part, "%Y%m%d")
    date_formatted = time.strftime("%d%b%Y", date_formatted)
    print("Entering loop for index {}: filename={}, date_part={}, date_formatted={}".format(index, filename, date_part, date_formatted))


    # Add the BandsExtractorOp node for the current slave file with the updated date, SWATH, and POLARISATION
    xml_content += '''
<node id="BandsExtractor_{}">
  <operator>BandsExtractorOp</operator>
  <sources>
    <sourceProduct refid="Read_{}"/>
  </sources>
  <parameters>
    <sourceBandNames>i_{}_{}_slv1_{},q_{}_{}_slv1_{}</sourceBandNames>
  </parameters>
</node>
    '''.format(index, index, IWs, POLARISATION, date_formatted, IWs, POLARISATION, date_formatted)


# Create BandMerge node
xml_content += '''
<node id="BandMerge">
    <operator>BandMerge</operator>
    <sources>
      <sourceProduct refid="BandsExtractor_1"/>
'''

# Add sources for slave BandsExtractorOp nodes
for index in range(1, len(slave_files)-1):
    xml_content += '      <sourceProduct.{} refid="BandsExtractor_{}"/>\n'.format(index, index+1)

xml_content += '''
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <sourceBands/>
      <geographicError>1.0E-5</geographicError>
    </parameters>
  </node>
'''

# Create BandMathsOp node to calculate avgIntensity
xml_content += '''
<node id="BandMaths">
  <operator>BandMaths</operator>
  <sources>
    <sourceProduct refid="BandMerge"/> 
  </sources>
  <parameters>
    <targetBands>
      <targetBand>
        <name>avgIntensity</name>
        <type>float32</type>
        <expression>
'''

# Generate the expression based on the number of slave images
expression = []

# Add expressions for slave images (indexes 1 and onwards)
for index, file in enumerate(slave_files[1:], start=1):
    # Extract the date part from the filename in the file path
    filename = os.path.basename(file)  # Get just the filename
    date_part = filename.split("_")[1]  # Assuming the date is the second part between underscores
    # Convert the date format (e.g., "20190709" to "09Jul2019")
    date_formatted = time.strptime(date_part, "%Y%m%d")
    date_formatted = time.strftime("%d%b%Y", date_formatted)

    expression.append('i_{}_{}_slv1_{} * i_{}_{}_slv1_{} + q_{}_{}_slv1_{} * q_{}_{}_slv1_{}'.format(
        IWs, POLARISATION, date_formatted, IWs, POLARISATION, date_formatted,
        IWs, POLARISATION, date_formatted, IWs, POLARISATION, date_formatted))

expression_str = ' + '.join(expression)
xml_content += '(' + expression_str + ') / ({}-1)'.format(Nstack)

# Close the XML content
xml_content += '''
        </expression>
      </targetBand>
    </targetBands>
  </parameters>
</node>
'''

# Add the Subset node
xml_content += '''
  <node id="Subset">
    <operator>Subset</operator>
    <sources>
      <sourceProduct refid="BandMaths"/>
    </sources>
    <parameters class="com.bc.ceres.binding.dom.XppDomElement">
      <sourceBands/>
      <region>0,0,24608,2844</region>
      <geoRegion>POLYGON</geoRegion>
      <subSamplingX>1</subSamplingX>
      <subSamplingY>1</subSamplingY>
      <fullSwath>false</fullSwath>
      <tiePointGridNames/>
      <copyMetadata>true</copyMetadata>
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
      <demResamplingMethod>DEM_RESAMPLING</demResamplingMethod>
      <imgResamplingMethod>DEM_RESAMPLING</imgResamplingMethod>
      <pixelSpacingInMeter>13.95755</pixelSpacingInMeter>
      <pixelSpacingInDegree>1.2538280493862425E-4</pixelSpacingInDegree>
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
      <sourceProduct refid="Subset"/>
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

xml_file_path = os.path.join(graphs_dir, 'SEN_slave_intensity.xml')
with open(xml_file_path, 'w') as xml_file:
    xml_file.write(xml_content)

print('Generated XML saved as SEN_slave_intensity.xml')


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
    filedata = filedata.replace('Nstack',str(Nstack))
    filedata = filedata.replace('IWs',IWs)
    filedata = filedata.replace('POLARISATION',POLARISATION)
    filedata = filedata.replace('DEMNAME', DEMNAME)
    filedata = filedata.replace('DEMFILE', DEMFILE)
    filedata = filedata.replace('DEM_RESAMPLING', DEMRESAMPLING)
    filedata = filedata.replace('OUTPUTFILE_INT_TC', outputname_int_tc)
    filedata = filedata.replace('OUTPUTFILE_INT', outputname_int)
    filedata = filedata.replace('OUTPUTINTFOLDER', intensityfolder)
    filedata = filedata.replace('POLYGON',polygon)

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

