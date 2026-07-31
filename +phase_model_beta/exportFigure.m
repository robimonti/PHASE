function exportFigure(fig,pathValue)
%EXPORTFIGURE Reliably render a report figure before Python embeds it.

folder = fileparts(pathValue);
if ~isfolder(folder), mkdir(folder); end
drawnow;
try
    exportgraphics(fig,pathValue,'Resolution',300);
catch
    print(fig,pathValue,'-dpng','-r300');
end
drawnow;
if ~isfile(pathValue)
    error('PHASE_Model_beta:figureExportFailed', ...
        'MATLAB did not create the report figure: %s',pathValue);
end
info = dir(pathValue);
if isempty(info) || info.bytes == 0
    error('PHASE_Model_beta:figureExportFailed', ...
        'MATLAB created an empty report figure: %s',pathValue);
end
end
