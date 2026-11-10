function validateCskResume(folder,bounds,masterProcessing,firstStep)
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
    token=regexp(fileread(file),'(?s)<geoRegion>\s*POLYGON\s*\(\((.*?)\)\)\s*</geoRegion>','tokens','once');
    if isempty(token),error('PHASE:CskResumeProvenanceMissing','Cannot read the applied CSK AOI from %s.',file);end
    numbers=regexp(token{1},'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?','match');
    values=str2double(numbers);
    if mod(numel(values),2)~=0 || numel(values)<8 || any(~isfinite(values))
        error('PHASE:CskResumeProvenanceMissing','Invalid applied CSK polygon in %s.',file);
    end
    xy=reshape(values,2,[]).';
    applied=[min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))];
    matches=all(abs(applied-bounds)<=1e-9);
    % Older releases wrote only 3 decimals. Compare the exact old rounding,
    % rather than treating every nearby new AOI as the same selection.
    legacy=all(abs(applied-round(applied,3))<=1e-12);
    matches=matches || (legacy && all(abs(applied-round(bounds,3))<=1e-9));
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
