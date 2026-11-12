"""Terrain-aware CSK radar crops using the installed SNAP Java model."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import shutil
import xml.etree.ElementTree as ET
from phase_subprocess import live_popen
from phase_dem import prepare_external_dem

ALGORITHM = 'terrain-aware-v1'
# One registration window around the ground AOI, also covering interpolation edges.
REGISTRATION_HALO = 128


def snap_runtime(gpt):
    configured = os.environ.get('PHASE_SNAP_HOME')
    roots = [Path(configured)] if configured else list(Path(gpt).resolve().parents)[:4]
    for root in roots:
        modules = [root/name/'modules' for name in ('snap','microwavetbx','s1tbx')]
        if not (root/'snap/modules').is_dir():
            continue
        jars = [str(jar) for folder in modules if folder.is_dir()
                for jar in folder.rglob('*.jar') if 'gdal' not in jar.name.lower()]
        candidates = [root/'.install4j/jre.bundle/Contents/Home/bin/java',
                      root/'jre/bin/java.exe', root/'jre/bin/java',
                      root/'jre.bundle/Contents/Home/bin/java']
        candidates += list(root.glob('.install4j/jre*/bin/java*'))
        if os.environ.get('JAVA_HOME'):
            candidates += [Path(os.environ['JAVA_HOME'])/'bin/java.exe',Path(os.environ['JAVA_HOME'])/'bin/java']
        if shutil.which('java'):candidates.append(Path(shutil.which('java')))
        java = next((path for path in candidates if path.is_file()), None)
        if java and jars:
            return str(java), jars
    raise RuntimeError('Cannot locate the SNAP bundled Java runtime and microwave modules. Select gpt from the complete SNAP installation.')



def classpath_jar(helper, jars, project):
    """Keep the Java command below Windows' limit; remove identical mirrored JARs."""
    import zipfile
    items=[Path(jar) for jar in jars]
    key=hashlib.sha256(json.dumps([(str(p),p.stat().st_size,p.stat().st_mtime_ns) for p in items]
                                 +[(str(helper),0,0)]).encode()).hexdigest()[:20]
    folder=Path(project)/'auxiliary'/'SNAP';folder.mkdir(parents=True,exist_ok=True)
    output=folder/f'classpath_{key}.jar'
    if output.is_file():return output
    seen=set();urls=[helper.resolve().as_uri()+'/']
    from collections import Counter
    counts=Counter(p.name for p in items)
    for path in items:
        # Distinct versions with the same filename must retain SNAP's ordering.
        identity=(path.name,hashlib.sha256(path.read_bytes()).hexdigest()) if counts[path.name]>1 else (str(path),'')
        if identity in seen:continue
        seen.add(identity);urls.append(path.resolve().as_uri())
    line='Class-Path: '+' '.join(urls)
    folded='Manifest-Version: 1.0\r\n'+line[:70]+'\r\n'
    line=line[70:]
    while line:
        folded+=' '+line[:69]+'\r\n';line=line[69:]
    folded+='\r\n'
    temporary=output.with_suffix('.jar.tmp')
    with zipfile.ZipFile(temporary,'w') as archive:
        archive.writestr('META-INF/MANIFEST.MF',folded.encode('ascii'))
    os.replace(temporary,output)
    return output


def calculate_region(source, config, destination):
    required = ['DEMNAME','DEMFILE','LONMIN','LONMAX','LATMIN','LATMAX','GPTBIN_PATH','PROJECTFOLDER']
    missing = [key for key in required if key not in config]
    if missing:
        raise ValueError('Missing terrain-aware CSK configuration: '+', '.join(missing))
    dem_file, no_data = prepare_external_dem(config['DEMNAME'],config['DEMFILE'],config['PROJECTFOLDER'])
    java, jars = snap_runtime(config['GPTBIN_PATH'])
    helper = Path(__file__).resolve().parent.parent/'java'
    class_file = helper/'PhaseCskGeometry.class'
    if not class_file.is_file():
        raise RuntimeError('The PHASE terrain-aware CSK geometry helper is missing. Repair/update PHASE.')
    saved=json.loads((helper/'PhaseCskGeometry.build.json').read_text(encoding='utf-8'))
    if (saved['classSHA256']!=hashlib.sha256(class_file.read_bytes()).hexdigest()
            or saved['sourceSHA256']!=hashlib.sha256((helper/'PhaseCskGeometry.java').read_bytes().replace(b'\r\n',b'\n')).hexdigest()):
        raise RuntimeError('The CSK geometry helper does not match its source/checksum manifest. Repair/update PHASE.')
    command = [java,'-Xmx2G','-cp',str(classpath_jar(helper,jars,config['PROJECTFOLDER'])),
               'PhaseCskGeometry',str(source),config['DEMNAME'],dem_file,no_data,
               *[str(config[key]) for key in ('LONMIN','LONMAX','LATMIN','LATMAX')],
               str(REGISTRATION_HALO),str(destination)]
    print('Projecting the requested CSK ground AOI into radar coordinates using SNAP DEM and orbit...',flush=True)
    process = live_popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    process.communicate()
    if process.returncode != 0 or not Path(destination).is_file():
        raise RuntimeError(f'Terrain-aware CSK projection failed (exit {process.returncode}). See the SNAP geometry diagnostics above.')
    result = json.loads(Path(destination).read_text(encoding='utf-8'))
    if result.get('algorithm') != ALGORITHM or len(result.get('region',[])) != 4:
        raise RuntimeError('Invalid terrain-aware CSK geometry result.')
    dem_hash=None
    if config['DEMNAME'].strip().lower()=='external dem':
        with Path(config['DEMFILE']).open('rb') as stream:
            dem_hash=hashlib.file_digest(stream,'sha256').hexdigest()
    result.update(demSourceSHA256=dem_hash,source=str(Path(source).resolve()),demFile=dem_file,demName=config['DEMNAME'],demSourceFile=config['DEMFILE'],demNoData=float(no_data),
                  helperSHA256=hashlib.sha256(class_file.read_bytes()).hexdigest())
    Path(destination).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def subset_product(source, output, config, graph_path):
    """Calculate and apply a pixel subset; commit provenance only on success."""
    graph_path = Path(graph_path)
    pending = graph_path.with_suffix('.geometry.pending.json')
    manifest = graph_path.with_suffix('.geometry.json')
    pending.unlink(missing_ok=True)
    try:
        result = calculate_region(source,config,pending)
        graph = ET.parse(Path(config['GRAPHSFOLDER'])/'CSK_AOI_subset.xml')
        graph.find(".//node[@id='Read']/parameters/file").text = str(source)
        graph.find(".//node[@id='Write']/parameters/file").text = str(output)
        graph.find('.//region').text = ','.join(str(int(v)) for v in result['region'])
        # geoRegion overrides region in SNAP; clear it deliberately.
        graph.find('.//geoRegion').text = None
        graph.write(graph_path,encoding='utf-8',xml_declaration=True)
        command = [config['GPTBIN_PATH'],str(graph_path),'-c',config['CACHE'],'-q',config['CPU']]
        print(f'CSK terrain-aware radar crop: {result["region"]}; registration halo: {REGISTRATION_HALO} pixels.',flush=True)
        process = live_popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        process.communicate()
        if process.returncode != 0:
            raise RuntimeError(f'SNAP CSK master/slave subsetting failed (exit {process.returncode}). See the SNAP output above.')
        if not Path(output).is_file():
            raise RuntimeError(f'SNAP returned success but did not create the CSK master/slave subset: {output}')
        result['subsetGraphSHA256']=hashlib.sha256(graph_path.read_bytes()).hexdigest()
        result['outputHeaderSHA256']=hashlib.sha256(Path(output).read_bytes()).hexdigest()
        result['output']=str(Path(output).resolve())
        pending.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        os.replace(pending,manifest)
        return result
    finally:
        pending.unlink(missing_ok=True)


def finalize_export(project, export, master):
    """Record the ground AOI and measure coverage from the actual SNAP export."""
    import numpy as np
    project,export=Path(project),Path(export)
    provenance=project/'graphs/CSK_subsetMgraph2run.geometry.json'
    if not provenance.is_file():
        raise RuntimeError('This CSK export has no terrain-aware crop provenance. Rebuild from Step 1 with master processing enabled.')
    geometry=json.loads(provenance.read_text(encoding='utf-8'))
    if geometry.get('algorithm')!=ALGORITHM:
        raise RuntimeError('Unsupported CSK crop provenance. Rebuild the preprocessing.')
    west,east,south,north=geometry['requestedBounds']
    lon=export/'geo'/f'{master}.lon';lat=export/'geo'/f'{master}.lat'
    if not lon.is_file() or not lat.is_file() or lon.stat().st_size!=lat.stat().st_size:
        raise RuntimeError('Missing or inconsistent SNAP geolocation export.')
    x=np.memmap(lon,dtype='>f4',mode='r');y=np.memmap(lat,dtype='>f4',mode='r')
    if 'region' in geometry and len(x)!=int(geometry['region'][2])*int(geometry['region'][3]):
        raise RuntimeError('SNAP export dimensions differ from the terrain-aware master crop. Rebuild the CSK pipeline before PSI.')
    valid_count=inside_count=0;occupied=np.zeros((32,32),dtype=bool)
    for start in range(0,len(x),500000):
        a=x[start:start+500000];b=y[start:start+500000]
        valid=np.isfinite(a)&np.isfinite(b)&(np.abs(a)<=180)&(np.abs(b)<90)&~((a==0)&(b==0))
        inside=valid&(a>=west)&(a<=east)&(b>=south)&(b<=north)
        valid_count+=int(valid.sum());inside_count+=int(inside.sum())
        cols=np.minimum(31,((a[inside]-west)/(east-west)*32).astype(int))
        rows=np.minimum(31,((b[inside]-south)/(north-south)*32).astype(int))
        occupied[rows,cols]=True
    if not inside_count:
        raise RuntimeError('The DEM-corrected CSK export has no valid pixel inside the requested ground AOI. No usable PSI input was produced.')
    report=dict(schema=1,algorithm=ALGORITHM,requestedBounds=geometry['requestedBounds'],
                totalPixels=int(len(x)),validPixels=valid_count,pixelsInsideAOI=inside_count,
                coveredGeographicCells=int(occupied.sum()),geographicCells=1024,
                coverageNote='Descriptive geographic occupancy; terrain layover, shadow and scene boundaries can leave gaps. Not a PS/coherence estimate.')
    target=export/'phase_preprocessing_aoi.json'
    temporary=target.with_suffix('.json.tmp');temporary.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    os.replace(temporary,target)
    print(f'CSK ground AOI verified: {inside_count:,} valid radar pixels; {occupied.sum()}/1024 geographic coverage cells. PSI preparation will exclude the registration halo.',flush=True)
    return report


def save_products(folder, paths):
    """Publish the exact successful products; leave unrelated old files untouched."""
    folder=Path(folder)
    records=[]
    for path in paths:
        path=Path(path)
        if not path.is_file():raise RuntimeError(f'SNAP did not create the expected CSK product: {path}')
        records.append(dict(file=path.name,headerSHA256=hashlib.sha256(path.read_bytes()).hexdigest()))
    if not records:raise RuntimeError(f'No CSK products were produced in {folder}.')
    target=folder/'phase_products.json';temporary=target.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(dict(schema=1,algorithm=ALGORITHM,products=records),indent=2)+'\n',encoding='utf-8')
    os.replace(temporary,target)


def product_paths(folder):
    """Never mix stale products from a different master/AOI with a new run."""
    folder=Path(folder);manifest=folder/'phase_products.json'
    if not manifest.is_file():raise RuntimeError(f'Missing terrain-aware CSK product inventory in {folder}. Rebuild preprocessing from Step 1 with master processing enabled.')
    saved=json.loads(manifest.read_text(encoding='utf-8'))
    if saved.get('schema')!=1 or saved.get('algorithm')!=ALGORITHM or not saved.get('products'):
        raise RuntimeError(f'Invalid CSK product inventory: {manifest}')
    result=[]
    for item in saved['products']:
        name=item['file']
        if Path(name).name!=name or not name.endswith('.dim'):raise RuntimeError('Unsafe CSK product inventory entry.')
        path=folder/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=item['headerSHA256']:
            raise RuntimeError(f'CSK product is missing or changed: {path}. Rebuild the preceding processing stage.')
        result.append(str(path))
    return result
