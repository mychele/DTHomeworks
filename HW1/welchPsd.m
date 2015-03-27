function [ P_welch ] = welchPsd( z1, window, S )
% WELCHPSD This function computes the Welch estimation of the PSD of a
% random process.
%
%   z1: signal for which to compute the estimate
%   window: an array containing the window to use to frame the signal
%   S: number of overlapping samples 

D = length(window);
K = length(z1); % signal length
M_w = 1/D * sum(window.^2); % Power of the window
N_s = floor((K-D)/(D-S) + 1); % number of subsequences
P_per_w = zeros(K, N_s);
for s = 0:(N_s-1)
  z_s = window .* z1( s*(D-S) + 1 : s*(D-S) + D ); % 1.495 with index + 1
  Z_s = fft(z_s, K);
  P_per_w(:, s + 1) = abs(Z_s).^2/(D*M_w);
end
P_welch = sum(P_per_w, 2)/N_s;

end

