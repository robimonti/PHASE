function result = cachePsiMapTiles(launcherDir,uiDir,requests)
%CACHEPSIMAPTILES Reuse PHASE's verified tile downloader for uihtml on macOS.
if isempty(requests),result=struct('keys',{{}},'failedCount',0);return;end
if numel(requests)>120,error('PHASE:MapTileRequestLimit','Too many map tiles requested.');end
if isstruct(requests),requests=num2cell(requests);end % Keep singleton JSON arrays.
input=[tempname '.json'];output=[tempname '.json'];
cleanup=onCleanup(@()removeFiles(input,output)); %#ok<NASGU>
fid=fopen(input,'w');if fid<0,error('PHASE:MapTileIO','Cannot create tile request.');end
fprintf(fid,'%s',jsonencode(struct('requests',{requests})));fclose(fid);
python=phase_preprocessing_beta.resolvePython(getenv('PHASE_PYTHON'));
args={python,fullfile(fileparts(launcherDir),'pythonScripts','cache_phase_map_tiles.py'), ...
    '--cache-root',fullfile(uiDir,'map_tiles'),'--requests',input,'--output',output};
if ispc
    args=cellfun(@(v)['"' strrep(v,'"','""') '"'],args,'UniformOutput',false);
else
    args=cellfun(@phase_stamps_beta.quotePosix,args,'UniformOutput',false);
end
[status,message]=system(strjoin(args,' '));
if status~=0||~isfile(output),error('PHASE:MapTileDownload','%s',message);end
data=jsondecode(fileread(output));keys={};
if ~isempty(data.successful),keys={data.successful.key};end
result=struct('keys',{keys},'failedCount',numel(data.failed));
end

function removeFiles(varargin)
for k=1:nargin,if isfile(varargin{k}),delete(varargin{k});end,end
end
