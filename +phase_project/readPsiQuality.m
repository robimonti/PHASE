function data = readPsiQuality(exportFolder)
%READPSIQUALITY Load current screening results; never silently reuse stale data.
data=struct('available',false,'message','Quality screening has not been generated for this dataset.');
if isempty(exportFolder),return;end
path=fullfile(exportFolder,'quality','quality.json');
if ~isfile(path),return;end
try
    report=jsondecode(fileread(path));
    if report.schema~=1 || ~strcmp(report.algorithm,'psi-screening-v1')
        data.message='This quality report uses an unsupported format. Rebuild it in Preprocessing.';return;
    end
    validateReport(report);
    if isfield(report,'stackFiles')
        stack=dir(fullfile(exportFolder,'rslc','*.rslc'));subfolder='rslc';
        if isempty(stack),stack=dir(fullfile(exportFolder,'slc','*.slc'));subfolder='slc';end
        current=sort(string(arrayfun(@(f)[subfolder '/' f.name],stack,'UniformOutput',false)));
        expected=sort(string(report.stackFiles));
        if ~isequal(current(:),expected(:))
            data.message='The acquisition set changed after this quality map was generated. Rebuild it in Preprocessing.';return;
        end
    end
    project=fileparts(exportFolder);
    for k=1:numel(report.sources)
        source=report.sources(k);file=source.path;
        if isfield(source,'relativePath')
            relative=strrep(strrep(source.relativePath,'/',filesep),'\',filesep);
            file=fullfile(project,relative);
        end
        found=dir(file);
        if isempty(found)
            if source.requiredForFreshness
                data.message='The stack required to verify this quality map is missing. Rebuild it in Preprocessing.';return;
            end
            continue; % Optional SNAP products may have been deliberately cleaned.
        end
        epoch=double(java.io.File(file).lastModified())/1000;
        if found(1).bytes~=source.bytes || abs(epoch-source.modified)>.01
            data.message='The preprocessing data changed after this quality map was generated. Rebuild it in Preprocessing.';return;
        end
    end
    data=report;data.available=true;data.folder=fullfile(exportFolder,'quality');
    data.message='Preliminary screening only. Favourable areas are not guaranteed to produce PS.';
catch ME
    data.message=['Quality report unavailable: ' ME.message];
end
end

function validateReport(report)
g=report.grid;
assert(isnumeric(g.bounds)&&numel(g.bounds)==4&&all(isfinite(g.bounds))&& ...
    g.bounds(1)<g.bounds(2)&&g.bounds(3)<g.bounds(4), ...
    'PHASE:QualityReportInvalid','Invalid quality map bounds.');
assert(isscalar(g.rows)&&isscalar(g.columns)&&g.rows>=1&&g.columns>=1&& ...
    g.rows<=96&&g.columns<=96&&fix(g.rows)==g.rows&&fix(g.columns)==g.columns, ...
    'PHASE:QualityReportInvalid','Invalid quality grid dimensions.');
for name={'suitability','geometry','stability','coverage','coherence'}
    item=report.layers.(name{1});v=item.values;
    assert(isnumeric(v)&&numel(v)==g.rows*g.columns, ...
        'PHASE:QualityReportInvalid','Quality layer dimensions do not match the grid.');
    finite=v(isfinite(v));
    if strcmp(name{1},'suitability'),ok=all(ismember(finite,0:3));else,ok=all(finite>=0&finite<=1);end
    assert(ok,'PHASE:QualityReportInvalid','Quality layer contains invalid values.');
end
assert(isnumeric(report.reasonCodes)&&numel(report.reasonCodes)==g.rows*g.columns&& ...
    all(isfinite(report.reasonCodes),'all')&&all(ismember(report.reasonCodes,0:63),'all'), ...
    'PHASE:QualityReportInvalid','Invalid quality reason codes.');
assert(isstruct(report.reasonBits)&&all(isfield(report.reasonBits,{'bit','message'})), ...
    'PHASE:QualityReportInvalid','Invalid quality reason descriptions.');
end
