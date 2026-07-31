function [resolved, info] = resolveGpt(configured)
%RESOLVEGPT Resolve a SNAP gpt executable or bin folder cross-platform.

configured = strtrim(char(string(configured)));
if isempty(configured)
    error('PHASE:SnapGptMissing','SNAP gpt path is empty.');
end

candidates = {configured};
[folder,name,extension] = fileparts(configured);
if isfolder(configured)
    candidates = [candidates, {fullfile(configured,'gpt'), ...
        fullfile(configured,'gpt.exe')}]; %#ok<AGROW>
elseif isempty(extension)
    candidates = [candidates, {[configured '.exe']}]; %#ok<AGROW>
elseif strcmpi(extension,'.exe') && ~isempty(folder)
    candidates = [candidates, {fullfile(folder,name)}]; %#ok<AGROW>
end

resolved = configured;
exists = false;
for k = 1:numel(candidates)
    candidate = candidates{k};
    if isfile(candidate)
        resolved = char(java.io.File(candidate).getCanonicalPath());
        exists = true;
        break
    end
end

info = struct('configured',configured,'resolved',resolved, ...
    'exists',exists,'changed',~strcmp(configured,resolved));
end
