function queryPointsSelfTest()
%QUERYPOINTSSELFTEST GUI/TXT equivalence and input validation.
folder=tempname;mkdir(folder);cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
c=phase_model_beta.defaultConfig();c.flag_tsExtr=true;c.queryPointSource='manual';
c.queryPoints='[{"id":"GNSS_01","longitude":152.20890917,"latitude":-4.23823009},{"id":"GNSS_02","longitude":152.21,"latitude":-4.24}]';
points=phase_model_beta.parseQueryPoints(c.queryPoints);
[~,saved]=phase_model_beta.saveConfig(folder,c);
fromFile=readtable(saved.filepath_EXTR,'TextType','string');
assert(isequal(string(fromFile.id),points.id));
assert(max(abs(fromFile.longitude-points.longitude))<1e-12);
assert(max(abs(fromFile.latitude-points.latitude))<1e-12);
[loaded,~]=phase_model_beta.loadConfig(folder);assert(phase_model_beta.configsEqual(saved,loaded));
phase_model_beta.saveConfig(folder,saved);assert(isfile([saved.filepath_EXTR '.bak']));
bad={ '[]', '[{"id":"a","longitude":181,"latitude":0}]', ...
    '[{"id":"a","longitude":0,"latitude":null}]', ...
    '[{"id":"../a","longitude":0,"latitude":0}]', ...
    '[{"id":"a","longitude":0,"latitude":0},{"id":"A","longitude":1,"latitude":1}]'};
for k=1:numel(bad)
    failed=false;try,phase_model_beta.parseQueryPoints(bad{k});catch,failed=true;end
    assert(failed,'Invalid query points accepted.');
end
disp('Manual query points passed: TXT round trip, persistence, duplicate/range/path validation.');
end
