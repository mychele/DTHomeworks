function [ a, sigma_w ] = arModel( N, autoc )
% ARMODEL of order N, given the unbiased/biased estimate of the
% autocorrelation of the signal whose PSD has to be estimated

row1 = conj(autoc);
% create the Toeplitz R matrix
R = toeplitz(row1(1:N));
% create r vector
r = autoc(2:N+1);
% Yule-Walker equations
a = -inv(R)*r;
sigma_w = abs(autoc(1) + r'*a);  % Abs to correct rounding errors

end

