function value = geoSplinterJobPath(folder,filename,jobFolder)
%GEOSPLINTERJOBPATH Keep native job entries short and independent of install path.
% geoSplinter's legacy parser is unreliable with long absolute paths. The
% four processing folders are siblings, and runGeoSplinter starts the native
% process in their common parent.

folder = char(string(folder));
jobFolder = char(string(jobFolder));
folderFile = java.io.File(folder);
jobFile = java.io.File(jobFolder);
if ~folderFile.isAbsolute(), folderFile = java.io.File(pwd,folder); end
if ~jobFile.isAbsolute(), jobFile = java.io.File(pwd,jobFolder); end
folderCanonical = char(folderFile.getCanonicalPath());
jobCanonical = char(jobFile.getCanonicalPath());
[folderParent,folderName] = fileparts(folderCanonical);
jobParent = fileparts(jobCanonical);
if ispc
    sameParent = strcmpi(folderParent,jobParent);
else
    sameParent = strcmp(folderParent,jobParent);
end
if ~sameParent || isempty(folderName)
    error('PHASE_Model_beta:geoSplinterFolderLayout', ...
        'geoSplinter data and job folders must share a parent: %s and %s', ...
        folderCanonical,jobCanonical);
end
value = fullfile('.',folderName,char(string(filename)));
end
