function [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, f, flagSymmetry, flagWhite)
% Function to convert from covariance function Cy to Power Spectral Density
% PSDy, properly accounting for the sampling rate.
%
% SYNTAX:
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, f)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, fs)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, f, flagSymmetry)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, fs, flagSymmetry)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, f, flagWhite)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(CPSDy, fs, flagWhite)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, f, flagSymmetry, flagWhite)
%      [Cy, tau, PSDy, f, fs, df] = psd2cov(PSDy, fs, flagSymmetry, flagWhite)
%
% INPUT:
% - PSDy         --> PSD for different frequencies [u/sqrt(Hz)]
% - f            --> frequency vector [Hz]
%                         or
% - fs           --> sampling frequency [Hz]
% - flagSymmetry --> '--autoSymmetrize' (default), '--Symmetrize', '--noSymmetrize'
%                    o '--autoSymmetrize' -> automatically understand if
%                                            the covariance has to be symmetrized
%                    o '--Symmetrize'     -> force symmetrization. This
%                                            will replace the given tau
%                    o '--noSymmetrize'   -> assume that the covariance is
%                                            correctly symmetrized with tau = 0 in
%                                            the "middle" (pos. N/2+1 or (N+1)/2)
% - flagWhite    --> '--autoWhite' (default), '--forceWhite', or '--forceColoured'
%                    o '--autoWhite'     -> automatically check if the noise 
%                                           is white or not (this could be 
%                                           approximated)
%                    o '--forceWhite'    -> white noise assumption
%                    o '--forceColoured' -> coloured noise assumption
%
% OUTPUT:
% - Cy    --> covariance function vector [u]
% - tau   --> time lag vector [s]
% - PSDy  --> rearranged Power Spectral Density [u / Hz]
% - f     --> rearranged frequency vector [Hz]
% - fs    --> used sampling frequency [Hz]
% - df    --> used frequency discretization step [Hz]
%
% REMARKS:
% - when bot covariance function and time lag vector are passed they must
%   be consistent. If not error will be raised.
% - the covariance function can be defined in three ways:
%   o on positive and negative time lags, with tau=0 in the "middle"
%   o on positive and negative time lags, with tau=0 as the first element
%     (meaning that after tau_max, there is -tau_max)
%   o on positive lags, with tau=0 as the first element. In this case the
%     covariance will be automatically symmetrized, as well as the time lag
%     vector
% - check for covariance symmetry are present. If they are not satisfied
%   the function rise an error
%
% See also: time2freq iddft
%
% v 1.0
% Lorenzo Rossi
% DICA - Politecnico di Milano
% 2022-09-26

% -------------------------------------
th = 1e-10;  % imaginary part threshold
th2 = 1e-14; % numerical accuracy threshold

% -------------------------------------
% input check
if nargin == 2
    flagWhite   = '--autoWhite';
    flagSymmetry = '--autoSymmetrize';
elseif nargin == 3
    if sum(strcmpi(flagSymmetry, {'--autoWhite', '--forceWhite', '--forceColoured'})) == 1
        flagWhite = flagSymmetry;
        flagSymmetry = '--autoSymmetrize';
    else
        flagWhite = '--autoWhite';
    end
elseif nargin == 4
    if sum(strcmpi(flagSymmetry, {'--autoWhite', '--forceWhite', '--forceColoured'})) == 1
        tmp = flagSymmetry;
        flagSymmetry = flagWhite;
        flagWhite = tmp;
        clear tmp
    end
end

if length(f) == 1
    fs = f;
    flagF = 0;
else
    df  = mean(diff(f));
    flagF = 1;
end

% verify if the time vector is row or column vector
flagColumn = 1;
if size(PSDy,1) < size(PSDy,2)
    PSDy = PSDy';
    flagColumn = 0;
end
if size(f,1) < size(f,2)
    f = f';
end

% number of data
N = length(PSDy);

% check even or odd number of data
if mod(N, 2) == 1
    isOdd  = true;
    isEven = false;
elseif mod(N, 2) == 0
    isOdd  = false;
    isEven = true;
end

% compute sampling frequency
if flagF == 1
    if isOdd
        fs = 2*(max(abs(f)) + df / 2);
    elseif isEven
        fs = 2*max(abs(f));
    end
end

% check for symmetry
isSym = 0;  % flag for detected symmetry
if strcmpi(flagSymmetry, '--autoSymmetrize')
    if isOdd
        % isSymC1 = sum(PSDy(1:N) ~= PSDy(N:-1:1)) == 0;     % f = 0 in the middle
        % isSymC2 = sum(PSDy(2:N) ~= PSDy(N:-1:2)) == 0;     % f = 0 in the origin
        isSymC1 = max(abs(PSDy(1:N) - PSDy(N:-1:1))./abs(PSDy(1:N))) <= th;     % f = 0 in the middle
        isSymC2 = max(abs(PSDy(2:N) - PSDy(N:-1:2))./abs(PSDy(2:N))) <= th;     % f = 0 in the origin
        if flagF == 1
            % isSymF1 = sum(f(1:N) ~= - f(N:-1:1)) == 0;      % f = 0 in the middle
            % isSymF2 = sum(f(2:N) ~= - f(N:-1:2)) == 0;      % f = 0 in the origin
            isSymF1 = max(abs(f(1:N) + f(N:-1:1))./abs(f(1:N))) <= th;      % f = 0 in the middle
            isSymF2 = max(abs(f(2:N) + f(N:-1:2))./abs(f(2:N))) <= th;      % f = 0 in the origin
        else
            isSymF1 = isSymC1;
            isSymF2 = isSymC2;
        end
        if isSymC1 && isSymF1
            isSym = 1;   % f = 0 in the middle (correct position)
        elseif isSymC2 && isSymF2
            isSym = 2;   % f = 0 in the origin (no need of fftshift)
        elseif (isSymC1 || isSymC2 || isSymF1 || isSymF2)
            error('Check consistency between frequency and PSD vectors!');
        end
        clear isSymC1 isSymC2 isSymT1 isSymT2
    elseif isEven
        % isSymC1 = sum(PSDy(2:N) ~= PSDy(N:-1:2)) == 0;        % f = 0 in N/2+1 (correct!) or f = 0 in the origin
        % isSymC2 = sum(PSDy(1:N-1) ~= PSDy(N-1:-1:1)) == 0;    % f = 0 in N/2 (not correctly centred)
        isSymC1 = max(abs(PSDy(2:N) - PSDy(N:-1:2))./abs(PSDy(2:N))) <= th;        % f = 0 in N/2+1 (correct!) or f = 0 in the origin
        isSymC2 = max(abs(PSDy(1:N-1) - PSDy(N-1:-1:1))./abs(PSDy(1:N-1))) <= th;    % f = 0 in N/2 (not correctly centred)
        if flagF == 1
            % isSymF1 = sum(f(2:N) ~= -f(N:-1:2)) == 0;         % f = 0 in the N/2+1 (correct!)
            % isSymF2 = sum(f(1:N-1) ~= -f(N-1:-1:1)) == 0;     % f = 0 in N/2 (not correctly centred)
            % isSymF3 = sum(f(2:N/2) ~= -f(N:-1:N/2+2)) == 0;   % f = 0 in the origin
            isSymF1 = max(abs(f(2:N) + f(N:-1:2))./abs(f(2:N))) <= th;           % f = 0 in the N/2+1 (correct!)
            isSymF2 = max(abs(f(1:N-1) + f(N-1:-1:1))./abs(f(1:N-1))) <= th;     % f = 0 in N/2 (not correctly centred)
            isSymF3 = max(abs(f(2:N/2) + f(N:-1:N/2+2))./abs(f(2:N/2))) <= th;   % f = 0 in the origin
        else
            % [~, idxC] = max(PSD);
            isSymF1 = isSymC1; % && (idxC == N/2+1);
            isSymF2 = isSymC2;
            isSymF3 = false; % isSymC1 && (idxC == 1);
            if isSymF1
                warning('Without the frequency vector the element corresponding to f=0 in the PSD cannot be detected. It is assumed to be in position N/2+1 and not in position 1!');
            end
        end
        % set the symmetry flag
        if isSymC1 && isSymF1
            isSym = 1;  % f = 0 in the middle (correct position)
        elseif isSymC2 && isSymF2
            isSym = 3;  % f = 0 in the middle (wrong position, shifted by one)
        elseif isSymC1 && isSymF3
            isSym = 2;  % f = 0 in the origin (no need of fftshift)
        elseif (isSymC1 || isSymC2 || isSymF1 || isSymF2 || isSymF3)
            error('Check consistency between time lag and covariance vectors!');            
        end
        clear isSymC1 isSymC2 isSymT1 isSymT2 isSymT3
    end
end

switch lower(flagSymmetry)
    case '--autosymmetrize'
        switch isSym
            case 0 % no symmetry detected
                warning('No symmetry detected: the covariance will be reflected and frequency vector ignored!');
                PSDy   = [flipud(PSDy(2:end)); PSDy];
                N      = length(PSDy);
                isOdd  = true;
                isEven = false;
                % flagF  = 1;
                % f = df * (-(N-1)/2 : (N-1)/2)';
            case 1 % f = 0 in the middle (correct position)
                % do nothing, it is fine!
            case 2
                PSDy = fftshift(PSDy);
            case 3
                PSDy = [PSDy(N); PSDy(1:N-1)];
        end
    
    case '--symmetrize'
        PSDy = [flipud(PSDy(2:end)); PSDy];
        N    = length(PSDy);
        isOdd = true;
        isEven = false;
        % flagF = 1;
    
    case '--nosymmetrize'
        % do nothing, it is fine!
end

% -------------------------------------
% re-build frequency vector (according to sorted PSDy vector)
[f, fs, df] = time2freq(1/fs, N);

% -------------------------------------
% check for white noise
if strcmpi(flagWhite, '--autoWhite')
    if sum(PSDy == PSDy(f==0)) == N
        flagWhite = '--forceWhite';
        warning('PSD is processed as white noise!');
    else
        flagWhite = '--forceColoured';
        warning('PSD is processed as coloured noise!');
    end
end

% -------------------------------------
% compute Covariance
switch lower(flagWhite)
    case '--forcewhite'
        Cy = fftshift(ifft(ifftshift(PSDy)));
    case '--forcecoloured'
        Cy = fftshift(idtft(PSDy, 1/fs));
end

if ~isreal(Cy)
    if sum(abs((abs(Cy)-abs(real(Cy)))./abs(real(Cy))) >= th & max(abs(real(Cy)))./abs(Cy) <= th2) == 0
        Cy = real(Cy);
        warning('Complex result: imaginary part of the covariance is lower than the threshold (%.1e) and has been neglected', th)
    else
        error('Something goes wrong: imaginary part of the covariance is larger than the threshold!');
    end
end

% compute tau vector
if isOdd
    tau = (-(N-1)/2:(N-1)/2)'*1/fs;
elseif isEven
    tau = (-(N)/2:(N-1)/2)'*1/fs;
end

% restore row vectors, if required
if flagColumn == 0
    Cy   = Cy';
    tau  = tau';
    PSDy = PSDy';
    f = f';
end
