function ensurePythonDependencies(installRoot,figure)
%ENSUREPYTHONDEPENDENCIES Finish upgrades applied by pre-dependency updaters.
% A managed release checks once per version/interpreter. Quality screening
% also checks its configured interpreter, which can differ between projects.
prefix=fileparts(installRoot);
if ~isfile(fullfile(prefix,'install.json')),return;end
progress=[];
showUI=~isempty(figure)&&isvalid(figure)&&strcmp(figure.Visible,'on');
try
    python=phase_preprocessing_beta.resolvePython(getenv('PHASE_PYTHON'));
    version=phase_project.displayVersion(installRoot);
    marker=fullfile(prefix,'python-dependencies.json');
    if isfile(marker)
        try
            previous=jsondecode(fileread(marker));
            if isfield(previous,'python')&&isfield(previous,'version')&& ...
                    strcmp(previous.python,python)&&strcmp(previous.version,version)
                return;
            end
        catch
            % An interrupted marker write must not suppress the real check.
        end
    end
    if showUI
        progress=uiprogressdlg(figure,'Title','PHASE update', ...
            'Message','Checking Python GeoTIFF support; installing missing packages...', ...
            'Indeterminate','on');
        drawnow;
    end
    script=fullfile(installRoot,'PHASE_Preprocessing','snap2stamps','bin','phase_python_dependencies.py');
    [status,output]=system(sprintf('"%s" "%s" --ensure 2>&1',python,script));
    if status~=0,error('PHASE:PythonDependencies','%s',strtrim(output));end
    fid=fopen(marker,'w');
    if fid<0,error('PHASE:DependencyMarker','Cannot save dependency check: %s',marker);end
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(struct('version',version,'python',python)));
catch ME
    warning('PHASE:PythonDependencies','Python dependency check failed: %s',ME.message);
    if ~isempty(progress)&&isvalid(progress),delete(progress);end
    if showUI
        uialert(figure,['GeoTIFF support could not be prepared. Check the internet connection ' ...
            'and Python permissions, then restart PHASE or retry Build quality map. ' ME.message], ...
            'Python dependencies unavailable');
    end
    return;
end
if ~isempty(progress)&&isvalid(progress),delete(progress);end
end
