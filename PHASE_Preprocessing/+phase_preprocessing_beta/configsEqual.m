function tf = configsEqual(left, right)
%CONFIGSEQUAL Compare normalized preprocessing configurations.

if isempty(left) || isempty(right), tf = false; return; end
left = phase_preprocessing_beta.configToUi(left);
right = phase_preprocessing_beta.configToUi(right);
tf = isequaln(orderfields(left), orderfields(right));
end
