"""Ground/radar crop orchestration, provenance and portable Java packaging."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys
import zipfile
import pytest

ROOT=Path(__file__).resolve().parents[1]
BIN=ROOT/'PHASE_Preprocessing/snap2stamps/bin'
sys.path.insert(0,str(BIN))
import phase_csk_geometry as geometry


def test_portable_helper_matches_its_source_and_manifest():
    folder=BIN.parent/'java';saved=json.loads((folder/'PhaseCskGeometry.build.json').read_text())
    assert hashlib.sha256((folder/'PhaseCskGeometry.java').read_bytes().replace(b'\r\n',b'\n')).hexdigest()==saved['sourceSHA256']
    raw=(folder/'PhaseCskGeometry.class').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==saved['classSHA256']
    assert raw[:4]==b'\xca\xfe\xba\xbe' and struct.unpack('>H',raw[6:8])[0]==52


def test_helper_source_compiles_without_snap_or_build_dependencies(tmp_path):
    compiler=shutil.which('javac')
    if not compiler:pytest.skip('Maintainer JDK unavailable')
    version=subprocess.run([compiler,'-version'],capture_output=True,text=True,check=True)
    flags=['-source','8','-target','8'] if '1.8.' in version.stdout+version.stderr else ['--release','8']
    subprocess.run([compiler,*flags,'-d',str(tmp_path),str(BIN.parent/'java/PhaseCskGeometry.java')],check=True,capture_output=True)
    assert (tmp_path/'PhaseCskGeometry.class').is_file()


def test_manifest_classpath_escapes_spaces_and_keeps_distinct_versions(tmp_path):
    helper=tmp_path/'helper with spaces';helper.mkdir()
    paths=[]
    for folder,value in [('a','same'),('b','same'),('c','different')]:
        directory=tmp_path/folder;directory.mkdir();p=directory/'library.jar';p.write_text(value);paths.append(str(p))
    target=geometry.classpath_jar(helper,paths,tmp_path/'project')
    with zipfile.ZipFile(target) as z:raw=z.read('META-INF/MANIFEST.MF')
    assert all(len(line)<=72 for line in raw.split(b'\r\n'))
    unfolded=raw.decode().replace('\r\n ','')
    assert '%20' in unfolded and Path(paths[0]).as_uri() in unfolded
    assert Path(paths[1]).as_uri() not in unfolded and Path(paths[2]).as_uri() in unfolded
    stamp=target.stat().st_mtime_ns
    assert geometry.classpath_jar(helper,paths,tmp_path/'project').stat().st_mtime_ns==stamp


@pytest.mark.parametrize('exit_code',[0,7])
def test_pixel_subset_clears_geographic_override_and_commits_only_success(tmp_path,monkeypatch,exit_code):
    source=tmp_path/'source & CSK.h5';source.write_bytes(b'fixture')
    output=tmp_path/'output.dim';graph=tmp_path/'subset.xml'
    config=dict(GRAPHSFOLDER=str(BIN.parent/'graphs'),GPTBIN_PATH='unused',CACHE='128M',CPU='1')
    def calculate(source,config,path):
        result=dict(algorithm=geometry.ALGORITHM,region=[100,200,300,400],requestedBounds=[8,9,45,46])
        path.write_text(json.dumps(result));return result
    monkeypatch.setattr(geometry,'calculate_region',calculate)
    class Process:
        returncode=exit_code
        def communicate(self):
            if not exit_code:output.write_text('<Dimap_Document/>')
    monkeypatch.setattr(geometry,'live_popen',lambda *a,**k:Process())
    if exit_code:
        with pytest.raises(RuntimeError,match='exit 7'):geometry.subset_product(source,output,config,graph)
        assert not graph.with_suffix('.geometry.json').exists()
    else:
        result=geometry.subset_product(source,output,config,graph)
        import xml.etree.ElementTree as ET
        parsed=ET.parse(graph)
        assert parsed.find('.//geoRegion').text is None
        assert parsed.find('.//region').text=='100,200,300,400'
        assert parsed.find(".//node[@id='Read']/parameters/file").text==str(source)
        assert result['outputHeaderSHA256']==hashlib.sha256(output.read_bytes()).hexdigest()
        assert graph.with_suffix('.geometry.json').is_file()
    assert not graph.with_suffix('.geometry.pending.json').exists()


def test_export_coverage_records_ground_aoi_and_refuses_empty_intersection(tmp_path):
    import numpy as np
    project=tmp_path/'processing';(project/'graphs').mkdir(parents=True)
    export=project/'INSAR_20200101';(export/'geo').mkdir(parents=True)
    proof=project/'graphs/CSK_subsetMgraph2run.geometry.json'
    proof.write_text(json.dumps(dict(algorithm=geometry.ALGORITHM,requestedBounds=[9.5,11.5,45.5,47.5])))
    lon=np.array([[9,10,11,12]]*3,dtype='>f4');lat=np.array([[45]*4,[46]*4,[47]*4],dtype='>f4')
    lon.tofile(export/'geo/20200101.lon');lat.tofile(export/'geo/20200101.lat')
    report=geometry.finalize_export(project,export,'20200101')
    assert report['pixelsInsideAOI']==4 and report['validPixels']==12
    assert json.loads((export/'phase_preprocessing_aoi.json').read_text())['requestedBounds']==[9.5,11.5,45.5,47.5]
    proof.write_text(json.dumps(dict(algorithm=geometry.ALGORITHM,requestedBounds=[70,71,10,11])))
    with pytest.raises(RuntimeError,match='no valid pixel inside'):geometry.finalize_export(project,export,'20200101')


def test_native_snap_projection_on_external_metadata(tmp_path):
    product=os.environ.get('PHASE_TEST_CSK_PRODUCT');dem=os.environ.get('PHASE_TEST_CSK_DEM')
    gpt=os.environ.get('PHASE_TEST_CSK_GPT');bounds=os.environ.get('PHASE_TEST_CSK_BOUNDS')
    if not all([product,dem,gpt,bounds]):pytest.skip('Optional real SNAP/CSK metadata fixture unavailable')
    west,east,south,north=json.loads(bounds)
    config=dict(PROJECTFOLDER=str(tmp_path),DEMNAME='External DEM',DEMFILE=dem,GPTBIN_PATH=gpt,
                LONMIN=west,LONMAX=east,LATMIN=south,LATMAX=north)
    report=geometry.calculate_region(product,config,tmp_path/'native.json')
    assert report['algorithm']==geometry.ALGORITHM and report['demSamples']>4
    x,y,w,h=report['region'];width,height=report['sourceSize']
    assert x>=0 and y>=0 and w>0 and h>0 and x+w<=width and y+h<=height
    assert report['insideSourceSamples']>0


def test_product_inventory_excludes_stale_files_and_detects_changed_headers(tmp_path):
    current=tmp_path/'current.dim';current.write_text('current')
    old=tmp_path/'old.dim';old.write_text('stale old AOI')
    geometry.save_products(tmp_path,[current])
    assert geometry.product_paths(tmp_path)==[str(current)]
    current.write_text('changed')
    with pytest.raises(RuntimeError,match='missing or changed'):geometry.product_paths(tmp_path)
    geometry.save_products(tmp_path,[current])
    current.unlink()
    with pytest.raises(RuntimeError,match='missing or changed'):geometry.product_paths(tmp_path)


def test_snap_runtime_discovers_windows_and_unix_bundled_java(tmp_path,monkeypatch):
    root=tmp_path/'SNAP install';modules=root/'snap/modules';modules.mkdir(parents=True)
    (modules/'snap-core.jar').touch();(modules/'snap-gdal.jar').touch()
    gpt=root/'bin/gpt.exe';gpt.parent.mkdir();gpt.touch()
    java=root/'jre/bin/java.exe';java.parent.mkdir(parents=True);java.touch()
    monkeypatch.delenv('PHASE_SNAP_HOME',raising=False)
    found,jars=geometry.snap_runtime(gpt)
    assert Path(found)==java and len(jars)==1 and Path(jars[0]).name=='snap-core.jar'
