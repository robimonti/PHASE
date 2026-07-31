classdef ProcessingEngine < phase_preprocessing_beta.LegacyEngine
    %PROCESSINGENGINE Beta-specific corrections around the extracted engine.

    properties
        GenerateCoherence = true
        GenerateLia = true
        ActiveProcess = []
    end

    methods
        function reprojectGeoTiffIfGeographic(app, filePath, epsgCode, bandIndex)
            % Reproject by inverse-mapping every target pixel and resampling
            % the source raster. The stable helper only changed the raster
            % reference, which did not geometrically warp the pixel values.
            if nargin < 4, bandIndex = []; end
            normalizedPath = lower(strrep(filePath,'\','/'));
            if (~app.GenerateCoherence && contains(normalizedPath,'/coherence/')) || ...
                    (~app.GenerateLia && contains(normalizedPath,'/lia/'))
                return
            end
            app.notifyBetaProgress(94,['Reprojecting ' shortFileName(filePath)]);
            phase_preprocessing_beta.reprojectGeoTiff(filePath,epsgCode,bandIndex);
            app.notifyBetaProgress(96,['Reprojected ' shortFileName(filePath)]);
        end

        function notifyBetaProgress(app, percentage, phase)
            if isempty(app.ExternalProgressCallback), return; end
            try
                app.ExternalProgressCallback(struct( ...
                    'percentage',double(percentage), ...
                    'phase',char(string(phase)), ...
                    'indeterminate',false));
            catch
            end
        end

        function stopped = forceStopActiveProcess(app)
            stopped = false;
            if isempty(app.ActiveProcess), return; end
            try
                if app.ActiveProcess.isAlive()
                    stopped = phase_preprocessing_beta.forceStopProcess( ...
                        app.ActiveProcess);
                end
            catch
            end
        end
    end
end

function name = shortFileName(pathValue)
[~,stem,extension] = fileparts(pathValue);
name = [stem extension];
end
