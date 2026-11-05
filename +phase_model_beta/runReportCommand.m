function [status,output]=runReportCommand(app,command)
%RUNREPORTCOMMAND Preserve all workbook data, but never report false success.
[status,output]=phase_model_beta.runCommandHidden(app,command,'Excel report formatting');
if status~=0
    app.ReportFormattingErrors{end+1}=sprintf('exit %d: %s',status,strtrim(output));
end
end
