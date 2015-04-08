function [ autoc ] = autocorrelation_biased( z1, N_corr )
%AUTOCORRELATION Summary of this function goes here
%   Detailed explanation goes here

K = length(z1);
autoc = zeros(N_corr + 1, 1);
% we should use the unbiased estimator pg 82 1.478
for n = 1:(N_corr + 1)
    d = z1(n:K);
    b = conj(z1(1:(K - n + 1)));
    c = K; % check this scaling factor
    autoc(n) = d.' * b / c;
end

end