function writeJson(pathValue, value)
%WRITEJSON Atomically write a small project metadata file.

folder = fileparts(pathValue);
if ~isfolder(folder)
    error('PHASE:MetadataFolderMissing','Folder does not exist: %s.',folder);
end
temporary = [tempname(folder) '.json'];
cleanup = onCleanup(@() removeTemporary(temporary)); %#ok<NASGU>
fid = fopen(temporary,'w','n','UTF-8');
if fid < 0
    error('PHASE:MetadataWriteFailed','Cannot create %s.',temporary);
end
try
    fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
    fclose(fid);
catch ME
    fclose(fid);
    rethrow(ME);
end
[ok,message] = movefile(temporary,pathValue,'f');
if ~ok
    error('PHASE:MetadataWriteFailed','Could not save %s: %s',pathValue,message);
end
end

function removeTemporary(pathValue)
if isfile(pathValue), delete(pathValue); end
end
