function themeLegacyEngine(engine)
%THEMELEGACYENGINE Give advanced native tools the PHASE visual language.

fig = engine.UIFigure;
fig.Name = 'PHASE · Preprocessing advanced workspace';
fig.Color = [1 1 1];
screen = get(groot, 'ScreenSize');
% The extracted advanced workspace uses the stable app's absolute layout.
% Keep its proven canvas size instead of stretching controls out of alignment.
width = min(1200, screen(3) - 80); height = min(675, screen(4) - 120);
fig.Position = [max(20,(screen(3)-width)/2), max(40,(screen(4)-height)/2), width, height];

if ispc, fontName = 'Segoe UI'; else, fontName = 'SF Pro Text'; end
components = findall(fig);
for k = 1:numel(components)
    component = components(k);
    try
        if isprop(component, 'FontName'), component.FontName = fontName; end
        if isa(component, 'matlab.ui.container.Panel') || isa(component, 'matlab.ui.container.Tab')
            component.BackgroundColor = [1 1 1];
        elseif isa(component, 'matlab.ui.control.Button')
            component.BackgroundColor = [0.94 0.96 1.00];
            component.FontColor = [53 101 207] / 255;
            component.FontWeight = 'bold';
        elseif isprop(component, 'BackgroundColor') && ...
                (isa(component, 'matlab.ui.control.EditField') || ...
                 isa(component, 'matlab.ui.control.NumericEditField') || ...
                 isa(component, 'matlab.ui.control.DropDown') || ...
                 isa(component, 'matlab.ui.control.TextArea'))
            component.BackgroundColor = [0.975 0.978 0.985];
        end
    catch
    end
end

primary = {'StartButton','StartButton_2','SearchASFButton', ...
    'DownloadSelectedButton','DownloadAllButton','DownloadRunUpdateImagesButton'};
for k = 1:numel(primary)
    try
        button = engine.(primary{k});
        button.BackgroundColor = [53 101 207] / 255;
        button.FontColor = [1 1 1];
    catch
    end
end
for name = {'StopButton','StopButton_2'}
    try
        button = engine.(name{1});
        button.BackgroundColor = [203 46 108] / 255;
        button.FontColor = [1 1 1];
    catch
    end
end
end
