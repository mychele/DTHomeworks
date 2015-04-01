close all;
%clear all;
clc;

%% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length
autoc = autocorrelation(z, length(z)/5);

%% High pass filter

z = filter(hp, 1, z);


%% Filter pg 137 BC
N = 36;
w0 = 2*pi*0.78; % 0.78 is probably our freq
vec =  w0 * (0:N-1);
E_w0 = exp(1i * vec).';

A = sqrt(30); % see plot of real autoc, rule of thumb
B=A;
lambda = A^2/autoc(1); % autoc(1) should be sigma^2 of the noise, if white

gain = B/A; 

D = 10; % delay

c_opt = gain * lambda * exp(- 1i * w0 * D) * E_w0 / (1 + N*lambda);

[H, w] = freqz(c_opt, 1, 'whole');

figure, plot(w/(2*pi), 10*log10(H))
axis([0, 1, -30, 10])
figure, plot(w/(2*pi), angle(H))

spectra = filter(c_opt, 1, z);


%% Compute different spectral analysis

% PERIODOGRAM pg 84
Z = fft(spectra);
periodogr = abs(Z).^2/K;

% compute WELCH estimator pg 85 
D = 200; % window size
window = kaiser(D, 5.65);
S = D/2; %common samples
P_welch = welchPsd(spectra, window, S);

% CORRELOGRAM
N_corr = ceil(K/5); % N_corr is the order of the autocorrelation estimate
window_correlogram = kaiser(2*N_corr + 1, 5.65); % window centered around N_corr
correlogram = correlogramPsd(spectra, window_correlogram, N_corr);

% AR model: order estimation
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
autoc = autocorrelation(spectra, N_corr);
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
    [~, sigma_w(N)] = arModel(N, autoc);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))

% the knee is apparently at N = 3
% compute the vector of coefficients a
N = 3;
[a, sigma_w] = arModel(N, autoc);
[H, omega] = freqz(1, [1; a], K, 'whole');


clear a autoc D fir_bs_1 N N_corr S upp_limit window window_correlogram



%% Plot PSD estimate

figure
plot((1:K)/K, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
hold on
plot((1:K)/K, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
hold on
plot((1:K)/K, 10*log10(periodogr), 'c:')
hold on
plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
hold off
axis([0, 1, -40, 10])
legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)', 'Location', 'SouthEast')
title('Spectral analysis')