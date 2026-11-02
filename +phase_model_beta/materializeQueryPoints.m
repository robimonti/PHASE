function config=materializeQueryPoints(config,rootDir)
%MATERIALIZEQUERYPOINTS Persist GUI points in the existing TXT input format.
if ~config.flag_tsExtr || ~strcmp(config.queryPointSource,'manual'),return;end
points=phase_model_beta.parseQueryPoints(config.queryPoints);
destination=fullfile(rootDir,'query_points_gui.txt');
temporary=[tempname(rootDir) '.txt'];
cleanup=onCleanup(@() removeTemporary(temporary)); %#ok<NASGU>
writetable(points,temporary,'Delimiter','\t','Encoding','UTF-8');
if isfile(destination),copyfile(destination,[destination '.bak'],'f');end
[ok,message]=movefile(temporary,destination,'f');
if ~ok,error('PHASE_Model_beta:queryPointsSave','Cannot save GUI query points: %s',message);end
config.filepath_EXTR=destination;
end
function removeTemporary(file)
if isfile(file),delete(file);end
end
