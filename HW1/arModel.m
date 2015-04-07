function [ a, sigma_w ] = arModel( N, autoc )
%ARMODEL Summary of this function goes here
%   Detailed explanation goes here

row1 = conj(autoc);
%row1(1) = conj(row1(1));
R = toeplitz(row1(1:N));
r = autoc(2:N+1);
a = -inv(R)*r;
sigma_w = abs(autoc(1) + r'*a);  % Abs to correct rounding errors


end

