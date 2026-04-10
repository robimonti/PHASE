function [f, fs, df] = time2freq(ts, N)
% Function to define the frequency axis, given a time vector or the sampling
% time with the overall number of epochs.
%
% SYNTAX:
%         [f, fs, df] = time2freq(ts, N)
%         [f, fs, df] = time2freq(t)
%
% INPUT:
% - t  --> time vector [s]
%             or
% - ts --> sampling rate [s]
% - N  --> number of observations
%
% OUTPUT:
% - f  --> frequency vector [Hz]
% - fs --> sampling frequency [Hz]
% - df --> frequency discretization step [Hz]
%
% REMARK: it is suggested to start from the time vector to avoid possible
% misalignment due to approximated sampling rate ts
%
% v 1.0
% Lorenzo Rossi
% DICA - Politecnico di Milano
% 2022-09-22

if nargin == 1
    N = length(ts);
    ts = mean(unique(diff(ts)));
end

fs = 1 / ts;
df = fs / N;

if (mod(N,2) == 1)						% odd number of observations
   f = df * (-(N-1)/2 : (N-1)/2)';
else									% even number of observations
   f = df * (-N/2 : N/2-1)';
end