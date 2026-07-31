function pathValue = saveConfig(rootDir,config)
%SAVECONFIG Atomically save the standalone Model configuration.

if nargin < 1 || isempty(rootDir)
    rootDir = phase_model_beta.projectRoot();
end
pathValue = fullfile(rootDir,'input_model.mat');
backup = [pathValue '.bak'];
temporary = [tempname(rootDir) '.mat'];
cleanup = onCleanup(@() removeTemporary(temporary)); %#ok<NASGU>
save(temporary,'config','-mat');
if isfile(pathValue)
    copyfile(pathValue,backup,'f');
end
[ok,message] = movefile(temporary,pathValue,'f');
if ~ok
    error('PHASE_Model_beta:saveFailed','Could not save %s: %s',pathValue,message);
end
end

function removeTemporary(pathValue)
if isfile(pathValue)
    try, delete(pathValue); catch, end
end
end
