function cskResumeSelfTest()
%CSKRESUMESELFTEST Terrain-aware provenance, changed AOI and legacy rejection.
folder=tempname;mkdir(fullfile(folder,'graphs'));cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
bounds=[9.77528 9.8105398 45.0586109 45.0775897];
old=[9.718 9.823 45.044 45.087];writeGraph('CSK_subsetMgraph2run.xml',old);writeGraph('CSK_subsetgraph2run.xml',old);
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskResumeAreaChanged');
phase_preprocessing_beta.validateCskResume(folder,old,1,3);
% A full rebuild is allowed despite old provenance.
phase_preprocessing_beta.validateCskResume(folder,bounds,0,1);
% Geometry metadata retains the exact requested AOI.
writeGraph('CSK_subsetMgraph2run.xml',bounds);writeGraph('CSK_subsetgraph2run.xml',bounds);
phase_preprocessing_beta.validateCskResume(folder,bounds,1,3);
writeGraph('CSK_subsetgraph2run.xml',old);
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskResumeAreaChanged');
delete(fullfile(folder,'graphs','CSK_subsetMgraph2run.xml'));
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskResumeProvenanceMissing');
% Check real byte hashes against the known SHA-256 of 'abc'.
known='ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad';
writeGraph('CSK_subsetMgraph2run.xml',bounds);writeGraph('CSK_subsetgraph2run.xml',bounds);
for subfolder={'master','subset'}
    mkdir(fullfile(folder,subfolder{1}));fid=fopen(fullfile(folder,subfolder{1},'product.dim'),'w');fwrite(fid,'abc');fclose(fid);
end
dem=fullfile(folder,'dem.tif');fid=fopen(dem,'w');fwrite(fid,'abc');fclose(fid);
for name={'CSK_subsetMgraph2run.geometry.json','CSK_subsetgraph2run.geometry.json'}
    proof=struct('algorithm','terrain-aware-v1','requestedBounds',bounds,'output','old/location/product.dim', ...
        'outputHeaderSHA256',known,'demName','External DEM','demSourceSHA256',known);
    fid=fopen(fullfile(folder,'graphs',name{1}),'w');fprintf(fid,'%s',jsonencode(proof));fclose(fid);
end
phase_preprocessing_beta.validateCskResume(folder,bounds,1,3,'External DEM',dem);
fid=fopen(dem,'w');fwrite(fid,'changed');fclose(fid);
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3,'External DEM',dem),'PHASE:CskResumeDemChanged');
writeGraph('CSK_subsetMgraph2run.xml',bounds);
delete(fullfile(folder,'graphs','CSK_subsetMgraph2run.geometry.json'));
mustFail(@()phase_preprocessing_beta.validateCskResume(folder,bounds,1,3),'PHASE:CskTerrainRebuildRequired');
fprintf('CSK resume provenance regressions passed.\n');
    function writeGraph(name,b)
        fid=fopen(fullfile(folder,'graphs',name),'w');
        fprintf(fid,'<graph><geoRegion>POLYGON ((%.10f %.10f,%.10f %.10f,%.10f %.10f,%.10f %.10f,%.10f %.10f))</geoRegion></graph>', ...
            b(1),b(3),b(2),b(3),b(2),b(4),b(1),b(4),b(1),b(3));fclose(fid);
        metadata=struct('algorithm','terrain-aware-v1','requestedBounds',b);
        fid=fopen(fullfile(folder,'graphs',replace(name,'.xml','.geometry.json')),'w');
        fprintf(fid,'%s',jsonencode(metadata));fclose(fid);
    end
end
function mustFail(action,id)
try,action();catch ME,assert(strcmp(ME.identifier,id),'Unexpected error: %s',ME.message);return;end
error('Expected failure: %s',id);
end
