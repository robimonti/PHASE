function root = findRoot(folder)
%FINDROOT Find an enclosing PHASE 7 project from a data subfolder.

root = '';
folder = char(string(folder));
if isempty(folder), return; end
candidate = char(java.io.File(folder).getCanonicalPath());
while true
    if isfile(fullfile(candidate,'phase-project.json'))
        phase_project.open(candidate);
        root = candidate;
        return
    end
    parent = fileparts(candidate);
    if strcmp(parent,candidate), return; end
    candidate = parent;
end
end
