close all;
clear all;
clc;

%% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length


%% Compute different spectral analysis

% PERIODOGRAM pg 84
Z = fft(z);
periodogr = abs(Z).^2/K;

% compute WELCH estimator pg 85 
D = 200; % window size
window = kaiser(D, 5.65);
S = D/2; %common samples
P_welch = welchPsd(z, window, S);

% CORRELOGRAM
N_corr = ceil(K/5); % N_corr is the order of the autocorrelation estimate
window_correlogram = kaiser(2*N_corr + 1, 5.65); % window centered around N_corr
correlogram = correlogramPsd(z, window_correlogram, N_corr);

% AR model: order estimation
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
autoc = autocorrelation(z, N_corr);
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
axis([0, 1, -10, 40])
legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)', 'Location', 'SouthEast')
title('Spectral analysis')



%% Split continuous and spectral lines

firbp = -firbs;
firbp((length(firbs)+1)/2) = firbp((length(firbs)+1)/2) + 1;    % Bandpass
z_cont = filter(firbs, 1, z);
z_lines = filter(firbp, 1, z);

% Skip the transients
z_cont  = z_cont (length(firbp):length(z_cont));
z_lines = z_lines(length(firbp):length(z_lines));
K = length(z_cont);

% Periodogram of the continuous spectrum
CZ = fft(z_cont, K);
period_continuous = abs(CZ).^2/K;
figure
subplot(2,1,1);
plot((1:K)/K, 10*log10(period_continuous))
title('Continuous spectrum');
axis([0, 1, -5, 40])

% Plot the periodogram of the spectral lines
LZ = fft(z_lines, K);
period_lines = abs(LZ).^2/K;
subplot(2,1,2);
plot((1:K)/K, 10*log10(period_lines))
title('Spectral lines');
axis([0, 1, -5, 40])


%% AR continuous

% AR model: order estimation
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
N_corr = K/5;
autoc_cont = autocorrelation(z_cont, N_corr);
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
    [~, sigma_w(N)] = arModel(N, autoc_cont);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))
title('Variance for the AR model of the continuous part');

% the knee is apparently at N = 3
% compute the vector of coefficients a
N = 3;
[a_cont, sigma_w] = arModel(N, autoc_cont);
[H_cont, omega] = freqz(1, [1; a_cont], K, 'whole');

figure, plot(omega/(2*pi), 10*log10(sigma_w*abs(H_cont).^2), 'Color', 'm', 'LineWidth', 1);
title('AR model of the continuous part');

figure, zplane(roots([1; a_cont]))
title('Location of poles and zeros for the AR model of continuous part');


%% AR lines

% AR model: order estimation
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
N_corr = K/5;
autoc_lines = autocorrelation(z_lines, N_corr);
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
    [~, sigma_w(N)] = arModel(N, autoc_lines);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))
title('Variance for the AR model of the spectral lines');

% the knee is apparently at N = 2
% compute the vector of coefficients a
N = 2;
[a_lines, sigma_w] = arModel(N, autoc_lines);
[H_lines, omega] = freqz(1, [1; a_lines], K, 'whole');

figure, plot(omega/(2*pi), 10*log10(sigma_w*abs(H_lines).^2), 'Color', 'm', 'LineWidth', 1);
title('AR model of the spectral lines');

figure, zplane(roots([1; a_lines]))
title('Location of poles and zeros for the AR model of the spectral lines');
