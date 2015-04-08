function [ autoc ] = autocorrelation_biased( z1, N_corr )
%AUTOCORRELATION biased estimator pg 83 of Benvenuto Cherubini

K = length(z1);
autoc = zeros(N_corr + 1, 1);
for n = 1:(N_corr + 1)
    d = z1(n:K);
    b = conj(z1(1:(K - n + 1)));
    c = K;
    autoc(n) = d.' * b / c;
end

end