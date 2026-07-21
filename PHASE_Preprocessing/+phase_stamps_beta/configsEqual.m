function equal = configsEqual(left, right)
%CONFIGSEQUAL Compare configuration structs independent of field ordering.

leftNames = sort(fieldnames(left));
rightNames = sort(fieldnames(right));
if ~isequal(leftNames, rightNames)
    equal = false;
    return
end
equal = true;
for k = 1:numel(leftNames)
    name = leftNames{k};
    if ~isequaln(left.(name), right.(name))
        equal = false;
        return
    end
end
end
