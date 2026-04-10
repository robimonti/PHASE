function [Fy, f] = dtft(y, ts, h)
% Compute the Discrete-time Fourier Transform (DTFT) with sampling time ts
% and time translation tau = h * ts
%
% The convention for the data sorting is that y(t=0) = y(1)
% The output will be sorted with Fy(f=0) = Fy(N/2+1)     if N is even
%                                Fy(f=0) = Fy((N-1)/2+1) if N is odd
% namely, it is like applying fftshift after the standard Matlab fft
% command
%
% SYNTAX:
%   [Fy, f] = dtft(y, ts)
%   [Fy, f] = dtft(y, ts, h)
%
% INPUT:
%   y   = signal in time
%   ts  = sampling time [s]
%   h   = ts steps of time shift
%
% OUTPUT:
%   Fy  = signal in frequency
%   f   = frequency axis of the transform [Hz]
%
% SEE ALSO idtft
%
% Reference:
% https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform
%
% v 1.0
% Lorenzo Rossi
% DICA - Politecnico di Milano
% 2022-09-27

% compute the DTFT, recalling that ts y(nT) = y(n)
Fy = fftshift(fft(y * ts));

if nargin > 2 || nargout > 1
    % compute frequency vector
    N = length(y);
    f = time2freq(ts, N);
    
    % apply time shift (if required)
    if nargin > 2
        Fy = Fy .* exp(2i * pi * f * h * ts);
    end
end
