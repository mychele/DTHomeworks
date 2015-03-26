function [ autoc_complete ] = autocorrelation_complete( z1, N_corr )
%AUTOCORRELATION_COMPLETE Summary of this function goes here
%   Detailed explanation goes here

autoc = autocorrelation(z1, N_corr);

K = length(z1);
% consider formulas as more close as possible to the book
% make autocorrelation simmetric
autoc_complete = zeros(K, 1);
autoc_complete(1:N_corr + 1) = autoc;
temp = flipud(conj(autoc));
% it's the same as puttig the conjugate, flipped, at the end of this
% vector, since fft considers a periodic repetition of the signal
autoc_complete((K - N_corr + 1):K) = temp(1:length(temp)-1);

end

