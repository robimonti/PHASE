function report_path = diagnose_PHASE_StaMPS()
%DIAGNOSE_PHASE_STAMPS Capture the Windows step-6 runtime configuration.
% Run this from the same ASC_/DES_ folder as PHASE_StaMPS.mlapp, then send
% PHASE_StaMPS_diagnostic.txt to the maintainer.  The function is read-only
% apart from replacing its own report file.

report_path = fullfile(pwd, 'PHASE_StaMPS_diagnostic.txt');
if exist(report_path, 'file') == 2
    delete(report_path);
end
diary(report_path);
cleanup_diary = onCleanup(@() diary('off'));

fprintf('PHASE StaMPS diagnostic\n');
fprintf('Generated: %s\n', datestr(now, 31));
fprintf('MATLAB: %s\n', version);
fprintf('Computer: %s\n', computer);
fprintf('pwd: %s\n\n', pwd);

fprintf('Working-directory files\n');
fprintf('  PHASE_StaMPS.mlapp: %d\n', exist(fullfile(pwd, 'PHASE_StaMPS.mlapp'), 'file') == 2);
fprintf('  input_StaMPS.mat:   %d\n', exist(fullfile(pwd, 'input_StaMPS.mat'), 'file') == 2);
fprintf('  PATCH_1:            %d\n', exist(fullfile(pwd, 'PATCH_1'), 'dir') == 7);
fprintf('  psver.mat:          %d\n', exist(fullfile(pwd, 'psver.mat'), 'file') == 2);
fprintf('  parms.mat:          %d\n\n', exist(fullfile(pwd, 'parms.mat'), 'file') == 2);

installation_folder = '';
train_flag = NaN;
subtr_tropo = '';
tropo_method = '';
stamps_first_step = '';
stamps_last_step = '';
input_file = fullfile(pwd, 'input_StaMPS.mat');
if exist(input_file, 'file') == 2
    input_cfg = load(input_file, 'installation_folder', 'train_flag', ...
        'subtr_tropo', 'tropo_method', 'stamps_first_step', 'stamps_last_step');
    if isfield(input_cfg, 'installation_folder')
        installation_folder = input_cfg.installation_folder;
    end
    if isfield(input_cfg, 'train_flag'), train_flag = input_cfg.train_flag; end
    if isfield(input_cfg, 'subtr_tropo'), subtr_tropo = input_cfg.subtr_tropo; end
    if isfield(input_cfg, 'tropo_method'), tropo_method = input_cfg.tropo_method; end
    if isfield(input_cfg, 'stamps_first_step'), stamps_first_step = input_cfg.stamps_first_step; end
    if isfield(input_cfg, 'stamps_last_step'), stamps_last_step = input_cfg.stamps_last_step; end
end
fprintf('Configured StaMPS root: %s\n', installation_folder);
fprintf('Configured root exists: %d\n', ~isempty(installation_folder) && isfolder(installation_folder));
fprintf('Saved TRAIN enabled: %d\n', isequal(train_flag, 0));
fprintf('Saved subtr_tropo: %s\n', char(string(subtr_tropo)));
fprintf('Saved tropo_method: %s\n', char(string(tropo_method)));
fprintf('Saved StaMPS steps: %s -> %s\n', ...
    char(string(stamps_first_step)), char(string(stamps_last_step)));

expected_snaphu = '';
if ~isempty(installation_folder)
    expected_snaphu = fullfile(installation_folder, 'external', 'snaphu', 'bin', 'snaphu.exe');
end
fprintf('Expected snaphu.exe: %s\n', expected_snaphu);
fprintf('Expected snaphu exists: %d\n\n', ~isempty(expected_snaphu) && exist(expected_snaphu, 'file') == 2);

fprintf('Resolved MATLAB implementations\n');
fprintf('ps_unwrap -all:\n');
disp(which('ps_unwrap', '-all'));
fprintf('uw_stat_costs -all:\n');
disp(which('uw_stat_costs', '-all'));
fprintf('aps_weather_model -all:\n');
disp(which('aps_weather_model', '-all'));

ps_unwrap_impl = which('ps_unwrap');
old_unwrap = false;
if ~isempty(ps_unwrap_impl) && exist(ps_unwrap_impl, 'file') == 2
    ps_unwrap_src = fileread(ps_unwrap_impl);
    old_unwrap = ~isempty(regexp(ps_unwrap_src, ...
        '(?m)^\s*\[ph_uw_some\]\s*=\s*uw_nosnaphu\s*\(', 'once'));
end
fprintf('Active ps_unwrap uses obsolete uw_nosnaphu call: %d\n\n', old_unwrap);

if ispc
    [where_status, where_output] = system('where snaphu');
else
    [where_status, where_output] = system('command -v snaphu');
end
fprintf('OS lookup status: %d\n', where_status);
fprintf('OS lookup output:\n%s\n', where_output);

fprintf('\nDiagnosis\n');
if isempty(installation_folder) || ~isfolder(installation_folder)
    fprintf('FAIL: input_StaMPS.mat points to a missing/empty StaMPS installation folder.\n');
elseif isempty(expected_snaphu) || exist(expected_snaphu, 'file') ~= 2
    fprintf('FAIL: snaphu.exe was not installed in the expected StaMPS folder.\n');
elseif where_status ~= 0
    fprintf('FAIL: snaphu.exe exists but its bin directory is not on the MATLAB process PATH.\n');
elseif old_unwrap
    fprintf('FAIL: MATLAB selected an obsolete ps_unwrap.m known to loop on Windows.\n');
else
    fprintf('PASS: the static Windows step-6 runtime checks passed.\n');
end

clear cleanup_diary
fprintf('Diagnostic report saved to: %s\n', report_path);
end
