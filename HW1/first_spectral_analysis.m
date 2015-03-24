close all;
clear all;
clc;

% load data
z = load('data for hw1.mat');
% make a column vector
z1 = z.z';
% remove continuos component
z1 = z1 - mean(z1);

K = length(z1);

% compute autocorrelation with an estimator
autoc = zeros(K/5, 1);

for n = 1:K/5
    d = z1(n:K);
    b = conj(z1(1:(K - n + 1)));
    c = K - n + 1;
    autoc(n) = d.' * b / c;
end

% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
upp_limit = K/5 -1;
for N = 1:upp_limit
    
    row1 = conj(autoc);
    row1(1) = conj(row1(1));
    R = toeplitz(row1(1:N));
    
    a = -inv(R)*autoc(2:N+1);
    sigma_w(N) = autoc(1) + autoc(2:N+1)'*a;
    
    %     figure
    %     [H, omega] = freqz(1, [1; a], 'whole');
    %   plot(omega, 10*log(sigma_w*abs(H)));
    
end

figure
plot(1:upp_limit, 10*log10(sigma_w))

% the knee is apparently at N = 3
% compute the vector of coefficients a
% check it out, it is slightly scaled than pwelch 
N = 3;
row1 = conj(autoc);
row1(1) = conj(row1(1));
R = toeplitz(row1(1:N));

a = -inv(R)*autoc(2:N+1);
sigma_w = autoc(1) + autoc(2:N+1)'*a;

% compare different spectral analysis
figure
subplot(1, 3, 1)
[H, omega] = freqz(1, [1; a], 'whole');
plot(omega, 10*log10(sigma_w*abs(H)));
title('ar(3)');
subplot(1, 3, 2)
[welch, w] = pwelch(z1);
plot(w, 10*log10(welch))
title('welch');
subplot(1, 3, 3)
[period, w] = periodogram(z1);
plot(w, 10*log10(period))
title('periodogram');

