function validateCskResume(folder,bounds,masterProcessing,firstStep,demName,demFile)
%VALIDATECSKRESUME Refuse to reuse CSK crops made for a different AOI.
% Run before config/graph generation so original provenance is preserved.
if masterProcessing==0 && firstStep>1
    error('PHASE:InconsistentCskResume','Reprocessing the CSK master requires first preprocessing step 1.');
end
graphs={};
if masterProcessing~=0,graphs{end+1}='CSK_subsetMgraph2run.xml';end
if firstStep>=3,graphs{end+1}='CSK_subsetgraph2run.xml';end
for k=1:numel(graphs)
    file=fullfile(folder,'graphs',graphs{k});
    if ~isfile(file)
        error('PHASE:CskResumeProvenanceMissing', ...
            'Cannot verify the AOI of existing CSK crops (%s missing). Start from step 1 with master processing enabled.',graphs{k});
    end
    geometryFile=replace(file,'.xml','.geometry.json');
    if ~isfile(geometryFile)
        error('PHASE:CskTerrainRebuildRequired', ...
            'These CSK crops predate terrain-aware AOI selection. Start from step 1 with master processing enabled to rebuild them using DEM and orbit.');
    end
    saved=jsondecode(fileread(geometryFile));
    if ~isfield(saved,'algorithm')||~strcmp(saved.algorithm,'terrain-aware-v1')|| ...
            ~isfield(saved,'requestedBounds')||numel(saved.requestedBounds)~=4
        error('PHASE:CskResumeProvenanceMissing','Invalid terrain-aware crop metadata: %s.',geometryFile);
    end
    applied=double(saved.requestedBounds(:).');
    matches=all(abs(applied-bounds)<=1e-9);
    if isfield(saved,'subsetGraphSHA256')&&~strcmp(saved.subsetGraphSHA256,sha256(file))
        error('PHASE:CskResumeProvenanceMissing','CSK subset graph changed after processing. Rebuild from step 1 with master processing enabled.');
    end
    output='';
    if isfield(saved,'output')
        name=regexp(saved.output,'[^/\\]+$','match','once');
        if startsWith(graphs{k},'CSK_subsetM'),subfolder='master';else,subfolder='subset';end
        output=fullfile(folder,subfolder,name);
    end
    if isfield(saved,'output')&&(~isfile(output)|| ...
            (isfield(saved,'outputHeaderSHA256')&&~strcmp(saved.outputHeaderSHA256,sha256(output))))
        error('PHASE:CskResumeProvenanceMissing','The saved CSK subset is missing or changed. Rebuild it before resuming.');
    end
    if nargin>=6
        if ~isfield(saved,'demName')||~strcmp(saved.demName,demName)
            error('PHASE:CskResumeDemChanged','The DEM used for CSK cropping changed. Rebuild from step 1 with master processing enabled.');
        end
        if strcmpi(strtrim(demName),'External DEM')&& ...
                (~isfile(demFile)||~isfield(saved,'demSourceSHA256')||~strcmp(saved.demSourceSHA256,sha256(demFile)))
            error('PHASE:CskResumeDemChanged','The external DEM used for CSK cropping changed or is missing. Rebuild the master and slaves from step 1.');
        end
    end
    if ~matches
        values=num2cell([applied bounds firstStep]);
        message=sprintf(['The existing CSK crops use AOI [%.8f %.8f %.8f %.8f], but the selected AOI is [%.8f %.8f %.8f %.8f] ', ...
             '(lon min/max, lat min/max). Step %d reuses those crops and cannot apply the changed AOI. ', ...
             'Start preprocessing from step 1 with master processing enabled, or restore the original AOI.'], ...
            values{:});
        error('PHASE:CskResumeAreaChanged','%s',message);
    end
end
end

function value=sha256(file)
fid=fopen(file,'r');if fid<0,error('PHASE:CskResumeProvenanceMissing','Cannot read %s.',file);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
digest=java.security.MessageDigest.getInstance('SHA-256');
while true
    bytes=fread(fid,1048576,'*uint8');if isempty(bytes),break;end
    digest.update(typecast(bytes,'int8'));
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
