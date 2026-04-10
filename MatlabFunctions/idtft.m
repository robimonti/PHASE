function [y, t] = idtft(Fy, ts, h)
% Compute the inverse Discrete-time Fourier Transform (DTFT) with sampling
% time ts and time translation tau = h * ts
%
% The convention for the data sorting is Fy(f=0) = Fy(N/2+1)     if N is even
%                                        Fy(f=0) = Fy((N-1)/2+1) if N is odd
% namely, it is like applying fftshift after the standard Matlab fft
% command
% The output will be given as y(t=0) = y(1)
%
% SYNTAX:
%   [y, tau] = idtft(Fy, ts)
%   [y, tau] = idtft(Fy, ts, h)
%
% INPUT:
%   Fy  = signal in frequency
%   ts  = sampling time [s]
%   h   = ts steps of the time shift
%
% OUTPUT:
%   y   = signal in time
%   t   = time lag axis [s]
%
% SEE ALSO dtft
%
% Reference:
% https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform
%
% v 1.0
% Lorenzo Rossi
% DICA - Politecnico di Milano
% 2022-09-27

% compute frequency vector
N = length(Fy);
f = time2freq(N, ts);
    
if nargin > 2    
    % apply time shift (if required)
    Fy = Fy .* exp(-2i * pi * f * h * ts);
end

% compute the DTFT
y = ifft(ifftshift(Fy ./ ts));

if nargout > 1
    % time2lag(ts, N);
    if size(Fy,2) == 1
        t = (0:N-1)'*ts;
    else
        t = (0:N-1)*ts;
    end
end
