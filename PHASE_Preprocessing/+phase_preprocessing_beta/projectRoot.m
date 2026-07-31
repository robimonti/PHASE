function root = projectRoot()
%PROJECTROOT Absolute path of the PHASE repository containing the beta.

packageFolder = fileparts(mfilename('fullpath'));
preprocessingFolder = fileparts(packageFolder);
root = fileparts(preprocessingFolder);
end
