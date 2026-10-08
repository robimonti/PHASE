function [project, p] = open(projectRoot)
%OPEN Validate and read a PHASE project without changing its contents.

p = phase_project.paths(projectRoot);
if ~isfile(p.manifest)
    error('PHASE:ProjectMissing','No phase-project.json in %s.',p.root);
end
try
    project = jsondecode(fileread(p.manifest));
catch ME
    error('PHASE:ProjectManifestInvalid', ...
        'Cannot read %s: %s',p.manifest,ME.message);
end
if ~isstruct(project) || ~isfield(project,'schemaVersion') || ...
        ~isequal(project.schemaVersion,1) || ~isfield(project,'id') || ...
        ~isfield(project,'name') || ~isfield(project,'layout') || ...
        ~strcmp(char(string(project.layout)),'phase-project-v1')
    error('PHASE:ProjectManifestInvalid', ...
        'Unsupported or incomplete PHASE project manifest: %s.',p.manifest);
end
end
