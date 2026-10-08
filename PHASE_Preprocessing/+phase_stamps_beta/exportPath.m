function pathValue = exportPath(projectRoot, masterDate)
%EXPORTPATH Locate SNAP's INSAR_<date> export in old or new layouts.

projectRoot = char(string(projectRoot));
masterDate = char(string(masterDate));
if isempty(regexp(masterDate,'^\d{8}$','once'))
    error('PHASE_StaMPS_beta:invalidMasterDate', ...
        'Master date must be YYYYMMDD.');
end
if isfile(fullfile(projectRoot,'phase-project.json'))
    [~,p] = phase_project.open(projectRoot);
    pathValue = fullfile(p.preprocessing,['INSAR_' masterDate]);
else
    pathValue = fullfile(projectRoot,'PHASE_Preprocessing', ...
        ['INSAR_' masterDate]);
end
end
