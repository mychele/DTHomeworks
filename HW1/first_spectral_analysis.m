close all;
clear all;
clc;

% load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
fir_bs_1 = fir_bs_1.fir_bs_1;
% make a column vector
z1 = z.z.';

z1 = z1 - mean(z1);

K = length(z1);


% compute different spectral analysis

% PERIODOGRAM pg 84
Z = fft(z1);
periodogram_1 = abs(Z).^2/K;

% compute WELCH estimator pg 85 
D = 200; % window size
window = kaiser(D, 5.65); % rectangular, consider using a hamming or smth else
M_w = 1/D * sum(window.^2);
S = D/2; %common samples
N_s = floor((K-D)/(D-S) + 1); % number of subsequences
P_per_w = zeros(K, N_s);
for s = 0:(N_s-1)
  z_s = window .* z1( s*(D-S) + 1 : s*(D-S) + D ); % 1.495 with index + 1
  Z_s = fft(z_s, K);
  P_per_w(:, s + 1) = abs(Z_s).^2/(D*M_w);
end
P_welch = sum(P_per_w, 2)/N_s;


% compute autocorrelation with an estimator
N_corr = ceil(K/5);

autoc = zeros(N_corr + 1, 1);
% we should use the unbiased estimator pg 82 1.478
for n = 1:(N_corr + 1)
    d = z1(n:K);
    b = conj(z1(1:(K - n + 1)));
    c = K - n + 1; % check this scaling factor
    autoc(n) = d.' * b / c;
end

% consider formulas as more close as possible to the book
% make autocorrelation simmetric
autoc_complete = zeros(K, 1);
autoc_complete(1:N_corr + 1) = autoc;
temp = flipud(conj(autoc));
% it's the same as puttig the conjugate, flipped, at the end of this
% vector, since fft considers a periodic repetition of the signal
autoc_complete((K - N_corr + 1):K) = temp(1:length(temp)-1);


% CORRELOGRAM
window_correlogram = bartlett(2*N_corr + 1); % window centered around N_corr
window_complete = zeros(K, 1);
window_complete(1:N_corr + 1) = window_correlogram(N_corr + 1 : 2*N_corr + 1);
window_complete(K - N_corr + 1 : K) = window_correlogram(1 : N_corr);

windowed_autoc = autoc_complete.*window_complete;
correlogram = fft(windowed_autoc);

% ar model
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
upp_limit = 60;
for N = 1:upp_limit
    
    row1 = conj(autoc);
    row1(1) = conj(row1(1));
    R = toeplitz(row1(1:N));
    
    a = -inv(R)*autoc(2:N+1);
    sigma_w(N) = autoc(1) + autoc(2:N+1)'*a;
    
      %figure
      %[H, omega] = freqz(1, [1; a], 'whole');
      %plot(omega, 10*log(sigma_w*abs(H)));
    
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

[H, omega] = freqz(1, [1; a], 'whole');

figure
plot((1:K)/K, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
hold on
plot((1:K)/K, 10*log10(correlogram), 'Color', 'b', 'LineWidth', 1)
hold on
plot((1:K)/K, 10*log10(periodogram_1), 'c:')
hold on
plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
hold off
axis([0, 1, -5, 40])
legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)', 'Location', 'SouthEast')
title('Spectral analysis')


continuousz1 = filter(fir_bs_1, 1, z1);
% PERIODOGRAM pg 84
CZ = fft(continuousz1);
period_continuous = abs(CZ).^2/K;
figure
plot((1:K)/K, 10*log10(period_continuous), 'c:')
axis([0, 1, -5, 40])


