function valueType = chooseExportValueType(choice, trainEffective, lastStep, outputKind)
%CHOOSEEXPORTVALUETYPE Select one scientifically explicit ps_plot series.
choice = char(string(choice));
outputKind = char(string(outputKind));
if strcmp(outputKind,'wrapped')
    valueType = '';
    return
end
if ~strcmp(outputKind,'unwrapped')
    error('PHASE_StaMPS:exportOutputKind','Unsupported phase output: %s',outputKind);
end
switch choice
    case 'standard'
        valueType = 'v-do';
    case 'corrected'
        if trainEffective
            valueType = 'v-dao';
        elseif lastStep >= 8
            valueType = 'v-dso';
        else
            error('PHASE_StaMPS:correctionUnavailable', ...
                'Corrected export requires available TRAIN or completed StaMPS Step 8.');
        end
    otherwise
        error('PHASE_StaMPS:exportChoice','Unsupported export choice: %s',choice);
end
end
