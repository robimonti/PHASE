function filename = create_safe_filename(base_name)

% create_safe_filename Creates a safe filename from a given base name by replacing invalid
% characters with underscores.

% Characters to replace
invalidChars = [' ', '/', '\', ':', '*', '?', '"', '<', '>', '|', '&', '#', '%'];

% Replace invalid characters with underscores
filename = base_name;
for i = 1:length(invalidChars)
filename = strrep(filename, invalidChars(i), '_');
end

% Remove leading/trailing underscores and multiple consecutive underscores
filename = regexprep(filename,'^_|_$',''); %trim leading and trailing
filename = regexprep(filename,'_+','_'); %replace multiple with one

% Limit filename length
maxFilenameLength = 200;
if length(filename) > maxFilenameLength
  filename = filename(1:maxFilenameLength);
end

% Ensure the filename is not empty
if isempty(filename)
  filename = 'untitled';
end
end