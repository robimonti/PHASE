function [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, z, t, flagMode, flagFig)
% Compute the empirical cross-PSD of given signals y(t) and z(t).
%
% SYNTAX:
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, t)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, t, flagMode)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, t, flagFig)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, t, flagMode, flagFig)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, z, t)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, z, t, flagMode)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, z, t, flagFig)
%     [PSDyz, fyz, Cyz, tauyz] = sigs2epsd(y, z, t, flagMode, flagFig)
%
% INPUT:
% - y        --> first signal in time
% - z        --> second signal in time (if not present PSD of y is computed)
% - t        --> time vector [s]
% - flagMode --> define the padding to apply to the input. Possible options:
%                 o '--zeropadding' (default): apply zero padding (suggested)
%                 o '--noPadding':             do not apply padding (circular covariance computation)
% - flagFig  --> '--noFig' (default) or '--showFig'
% 
% OUTPUT:
% - PSDyz --> cross-PSD between y and z (or PSD of y, if z is not passed as input)
% - fyz   --> frequency axis [Hz]
% - Cyz   --> cross-covariance between y and z (or auto-covariance of y, if z is not passed as input)
% - tauyz --> time lag axis [s]
%
% See also: time2freq time2lag ddft psd2cov
%
% v 1.1
% Lorenzo Rossi
% DICA - Politecnico di Milano
% 2023-02-02
%
% Change log:
%
% v 1.1
% -----
% - psd2cov must be called with the '--noSymmetrize' flag to avoid probelm
%   in case of cross-PSDs

% input parser to use all the possible sintax
switch nargin
    case 2
        isCrossPSD = false;
        t = z;
        clear z
        flagMode = '--zeropadding';
        flagFig  = '--noFig';
        
    case 3
        if ischar(t)
            isCrossPSD = false;
            if sum(strcmpi(t, {'--zeropadding', '--nopadding'})) == 1
                flagMode = t;
                flagFig  = '--noFig';
            else
                flagFig = t;
                flagMode = '--zeropadding';
            end
            t = z;
            clear z
        else
           isCrossPSD = true;
           flagMode = '--zeropadding';
           flagFig  = '--noFig'; 
        end
            
    case 4
        if ischar(t) && ischar(flagMode)
            isCrossPSD = false;
            if sum(strcmpi(t, {'--zeropadding', '--nopadding'})) == 1
                flagFig  = flagMode;
                flagMode = t;
            elseif sum(strcmpi(t, {'--noFig', '--showFig'})) == 1
                flagFig = t;
            else
                error('Revise flags: they does not correspond to any possible option!');
            end
            t = z;
            clear z
        else
            isCrossPSD = true;
            if sum(strcmpi(flagMode, {'--zeropadding', '--nopadding'})) == 1
                flagFig  = '--noFig';
             elseif sum(strcmpi(flagMode, {'--noFig', '--showFig'})) == 1
                flagFig = flagMode;
                flagMode = '--zeropadding';
            else
                error('Revise flags: they does not correspond to any possible option!');
            end
        end
    case 5
        isCrossPSD = true;
        if sum(strcmpi(flagFig, {'--zeropadding', '--nopadding'})) == 1 && ...
                sum(strcmpi(flagMode, {'--noFig', '--showFig'})) == 1
            tmp = flagFig;
            flagFig  = flagMode;
            flagMode = tmp;
        elseif sum(strcmpi(flagMode, {'--zeropadding', '--nopadding'})) == 1 && ...
                sum(strcmpi(flagFig, {'--noFig', '--showFig'})) == 1
        else
            error('Revise flags: they does not correspond to any possible option!');
        end 
end

% number of elements
N = length(y);

% understanding the meaning of t
if length(t) == 1
    ts = t;
else
    ts  = mean(diff(t));
end
% fs = 1/ts;
% df = fs / N;

switch lower(flagMode)
    case '--zeropadding'
        y_pad  = [y; zeros(N-1, 1)];          % zero padding of y
        Fy = dtft(y_pad, ts);                 % Fourier transform of y

        if isCrossPSD
            z_pad  = [z; zeros(N-1, 1)];                % zero padding of z
            Fz = dtft(z_pad, ts);                       % Fourier transform of y
            PSDyz = (Fy .* conj(Fz)) ./ (N * ts);       % empirical cross-PSD
        else
            PSDyz = abs(Fy).^2 ./ (N * ts);             % empirical auto-PSD
        end
        fyz   = time2freq(ts, 2*N-1);                   % frequency axis
        
    case '--nopadding'
        Fy = dtft(y, ts);                         % Fourier transform of y
        if isCrossPSD
            Fz = dtft(z, ts);                     % Fourier transform of y
            PSDyz = (Fy .* conj(Fz)) ./ (N * ts); % empirical PSD
        else
            PSDyz = abs(Fy).^2 ./ (N * ts);       % empirical auto-PSD
        end            
        fyz   = time2freq(ts, N);
end

% compute covariance function in time
if nargout > 2 || strcmpi(flagFig, '--showFig')
    % [Cyz, tauyz] = psd2cov(PSDyz, fyz);
    [Cyz, tauyz] = psd2cov(PSDyz, fyz, '--noSymmetrize');
end

% figures
if strcmpi(flagFig, '--showFig')
    warning off
    figure;
    subplot(121)
    loglog(fyz, PSDyz, '.-');
    xlabel('Frequency [Hz]');
    ylabel('PSD_{yz}');
    title('Empirical PSD');
    
    subplot(122)
    plot(tauyz, Cyz, '.-');
    xlabel('Time [s]');
    ylabel('C_{yz}');
    title('Biased empirical covariance');
    warning on
end