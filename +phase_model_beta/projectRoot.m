function root = projectRoot()
%PROJECTROOT Absolute PHASE folder containing the standalone Model beta.

packageFolder = fileparts(mfilename('fullpath'));
root = fileparts(packageFolder);
end
