function folder=figureFolder(dataFolder)
%FIGUREFOLDER Route run diagnostics to the shared output figures directory.
[parent,name]=fileparts(dataFolder);
[runRoot,parentName]=fileparts(parent);
if strcmp(parentName,'files') && ismember(name,{'customModel','residualDiagnostics'})
    folder=fullfile(runRoot,'figures');
else
    folder=fullfile(dataFolder,'figures');
end
if ~isfolder(folder),mkdir(folder);end
end
