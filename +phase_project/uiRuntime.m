function uiDir = uiRuntime(sourceDir, module, cacheRoot)
%UIRUNTIME Stage UI files in a writable per-user cache.
% The installation remains read-only; generated map tiles and live logs live
% beside the staged HTML so uihtml can load them through relative URLs.

sourceDir = char(string(sourceDir));
module = char(string(module));
if ~isfolder(sourceDir) || isempty(regexp(module,'^[a-z0-9_-]+$','once'))
    error('PHASE:UiRuntimeInvalid','Invalid PHASE UI source or module name.');
end
if nargin < 3 || isempty(cacheRoot)
    cacheRoot = fullfile(prefdir,'PHASE','ui');
end
uiDir = fullfile(char(string(cacheRoot)),module);
ensureFolder(uiDir);
sourceFiles = dir(sourceDir);
for k = 1:numel(sourceFiles)
    name = sourceFiles(k).name;
    if sourceFiles(k).isdir || startsWith(name,'.'), continue; end
    [ok,message] = copyfile(fullfile(sourceDir,name),fullfile(uiDir,name),'f');
    if ~ok
        error('PHASE:UiRuntimeCopyFailed','Could not stage %s: %s',name,message);
    end
end
sourceAssets = fullfile(sourceDir,'assets');
if isfolder(sourceAssets)
    assetsDir = fullfile(uiDir,'assets');
    ensureFolder(assetsDir);
    assets = dir(sourceAssets);
    for k = 1:numel(assets)
        if assets(k).isdir || startsWith(assets(k).name,'.'), continue; end
        [ok,message] = copyfile(fullfile(sourceAssets,assets(k).name), ...
            fullfile(assetsDir,assets(k).name),'f');
        if ~ok
            error('PHASE:UiRuntimeCopyFailed','Could not stage asset %s: %s', ...
                assets(k).name,message);
        end
    end
end
ensureFolder(fullfile(uiDir,'runtime_logs'));
ensureFolder(fullfile(uiDir,'map_tiles'));
if ~isfile(fullfile(uiDir,'index.html'))
    error('PHASE:UiRuntimeMissingHtml','No index.html in %s.',uiDir);
end
end

function ensureFolder(pathValue)
if isfolder(pathValue), return; end
[ok,message] = mkdir(pathValue);
if ~ok, error('PHASE:UiRuntimeCreateFailed','%s',message); end
end
