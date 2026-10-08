function folder = stampsFolder(projectRoot)
%STAMPSFOLDER Parent directory for generated ASC/DES datasets.

projectRoot = char(string(projectRoot));
if isfile(fullfile(projectRoot, 'phase-project.json'))
    [~, paths] = phase_project.open(projectRoot);
    folder = paths.stamps;
    if ~isfolder(folder), mkdir(folder); end
else
    folder = projectRoot;
end
end
