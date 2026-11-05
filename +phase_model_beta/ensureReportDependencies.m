function ensureReportDependencies(app, executable)
%ENSUREREPORTDEPENDENCIES Fail before modelling if images cannot enter Excel.
script=fullfile(phase_model_beta.projectRoot(),'pythonScripts','reportDependencies.py');
command=sprintf('"%s" "%s" --ensure',executable,script);
[status,output]=phase_model_beta.runCommandHidden(app,command,'Excel report dependency check');
if status~=0
    error('PHASE_Model_beta:reportDependencies', ...
        'Excel report dependencies could not be prepared in %s. No modelling has started. Details: %s',executable,output);
end
end
