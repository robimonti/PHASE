function folder = dataFolder(rootDir)
%DATAFOLDER Resolve preprocessing data in PHASE 6 and PHASE 7 layouts.

rootDir = char(string(rootDir));
if isfile(fullfile(rootDir,'phase-project.json'))
    [~,p] = phase_project.open(rootDir);
    folder = p.preprocessing;
else
    folder = fullfile(rootDir,'PHASE_Preprocessing');
end
end
