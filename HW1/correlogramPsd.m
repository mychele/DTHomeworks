function [ correlogram ] = correlogramPsd( z1, window, N_corr )
% Compute the PSD estimate using the correlogram method

K = length(z1);
autoc_complete = autocorrelation_complete(z1, N_corr);
window_complete = zeros(K, 1);
window_complete(1:N_corr + 1) = window(N_corr + 1 : 2*N_corr + 1);
window_complete(K - N_corr + 1 : K) = window(1 : N_corr);

windowed_autoc = autoc_complete .* window_complete;
correlogram = fft(windowed_autoc);

end

