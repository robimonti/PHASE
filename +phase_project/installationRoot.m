function root = installationRoot()
%INSTALLATIONROOT Locate the installed PHASE code, independent of a project.

root = fileparts(fileparts(mfilename('fullpath')));
root = char(java.io.File(root).getCanonicalPath());
end
