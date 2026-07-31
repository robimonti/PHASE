function themeLegacyEngine(app)
%THEMELEGACYENGINE Apply the PHASE light visual language to the text engine.
% The scientific controls and callbacks remain those of the validated model.

phaseBlue = [53 101 207] ./ 255;
phasePink = [203 46 108] ./ 255;
graphite = [69 70 70] ./ 255;
surface = [0.985 0.988 0.995];

try
    app.UIFigure.Name = 'PHASE · Geospatial Model Beta';
    app.UIFigure.Color = surface;
catch
end

controls = findall(app.UIFigure);
for k = 1:numel(controls)
    control = controls(k);
    try
        if isprop(control,'FontName'), control.FontName = 'Helvetica'; end
        if isprop(control,'FontColor'), control.FontColor = graphite; end
        if isprop(control,'BackgroundColor')
            typeName = class(control);
            if contains(typeName,'Button')
                control.BackgroundColor = [0.94 0.96 1.0];
            elseif contains(typeName,'Panel') || contains(typeName,'Tab')
                control.BackgroundColor = [1 1 1];
            end
        end
    catch
    end
end

try
    app.StartButton.BackgroundColor = phaseBlue;
    app.StartButton.FontColor = [1 1 1];
catch
end
try
    app.SaveButton.BackgroundColor = [0.94 0.96 1.0];
    app.SaveButton.FontColor = phaseBlue;
catch
end
try
    app.LoadButton.BackgroundColor = [1.0 0.94 0.97];
    app.LoadButton.FontColor = phasePink;
catch
end
end
