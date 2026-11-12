function [bounds,signature]=preprocessingArea(exportRoot)
%PREPROCESSINGAREA Exact ground AOI associated with a terrain-aware CSK export.
bounds=[];signature=[];
file=fullfile(exportRoot,'phase_preprocessing_aoi.json');
if ~isfile(file)
    input=fullfile(fileparts(exportRoot),'input_preprocessing.mat');
    if isfile(input)
        vars=whos('-file',input);
        if ismember('constellation',{vars.name})
            cfg=load(input,'constellation');
            if strcmpi(char(string(cfg.constellation)),'CSK')
                error('PHASE:CskTerrainRebuildRequired', ...
                    'This CSK export predates terrain-aware AOI cropping. Rebuild preprocessing from step 1 with master processing enabled before preparing PSI again. Existing results are unchanged.');
            end
        end
    end
    return
end
saved=jsondecode(fileread(file));
if ~isfield(saved,'algorithm')||~strcmp(saved.algorithm,'terrain-aware-v1')|| ...
        ~isfield(saved,'requestedBounds')||numel(saved.requestedBounds)~=4
    error('PHASE:InvalidPreprocessingArea','Invalid CSK ground-area metadata: %s.',file);
end
bounds=double(saved.requestedBounds(:).');
if any(~isfinite(bounds))||bounds(1)>=bounds(2)||bounds(3)>=bounds(4)
    error('PHASE:InvalidPreprocessingArea','Invalid CSK ground-area bounds in %s.',file);
end
signature=struct('algorithm',saved.algorithm,'lonMin',bounds(1),'lonMax',bounds(2),'latMin',bounds(3),'latMax',bounds(4));
end
