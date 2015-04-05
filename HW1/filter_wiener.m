close all;
%clear all;
clc;

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
hp18 = load('hp18.mat');
hp18 = hp18.hp18.'; % make a column vector

z = z - mean(z); % remove average
K = length(z); % signal length
autoc = autocorrelation(z, length(z)/5);

%% High pass filter

zhp = filter(hp18, 1, z);

%% Filter pg 137 BC
N_wien = 31;             % Number of coefficients of the Wiener filter
w0 = 2*pi*0.78;     % 0.78 is probably our freq
E_w0 = exp(1i * w0 * (0 : N_wien-1)).';
A = sqrt(30); % see plot of real autoc, rule of thumb
B=A;
lambda = A^2/autoc(1); % autoc(1) should be sigma^2 of the noise, if white
gain = B/A;
D = 100;    % delay

c_opt = gain * lambda * exp(- 1i * w0 * D) * E_w0 / (1 + N_wien*lambda);


% --- Plot frequency response of Wiener filter

DTFTplot(c_opt, 50000);
title('Freq resp of Wiener filter')
ylim([-20 0])
% [H, w] = freqz(c_opt, 1, 'whole');
% figure, plot(w/(2*pi), 20*log10(abs(H)))
% axis([0, 1, -40, 5])
% figure, plot(w/(2*pi), angle(H))


% --- Filter using the Wiener filter
spectral_lines = filter(c_opt, 1, zhp);


% --- Plot original and filtered signal
% figure
% subplot(2,2,1)
% plot(real(z))
% title('Real part of original signal');
% subplot(2,2,3)
% plot(imag(z), 'r')
% title('Imaginary part of original signal');
% subplot(2,2,2)
% plot(real(spectra))
% title('Real part of spectral line');
% subplot(2,2,4)
% plot(imag(spectra), 'r')
% title('Imaginary part of spectral line');


% --- Plot spectrum of the spectral lines
plot_spectrum(spectral_lines, 1); % 1 is the order of the desired AR model
title('Spectral analysis of the signal after Wiener'), legend('Location', 'NorthWest')



% --- Recompute everything for AR

% Find the knee of sigma_w
%N_corr = length(spectral_lines)/5;
%autoc = autocorrelation(spectral_lines, N_corr);
%upp_limit = 60;
%sigma_w = zeros(1, upp_limit);
%for N = 1:upp_limit
%    [~, sigma_w(N)] = arModel(N, autoc);
%end
%figure, plot(1:upp_limit, 10*log10(sigma_w))  % Uncomment to plot sigma

% Choose order for AR and compute the vector of coefficients a
%N = 1;
%[a, sigma_w] = arModel(N, autoc);
%[H, omega] = freqz(1, [1; a], K, 'whole');
%figure, plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
%axis([0, 1, -40, 10])



% --- Continuous part of the signal

% Compute the opposite of the Wiener filter we just used
neg_wiener = -c_opt;
neg_wiener((N_wien-1)/2) = neg_wiener((N_wien-1)/2) + 1;  % The order is N_wien-1

% Filter original signal (not the HP'ed one!)
continuous = filter(neg_wiener, 1, z);

% Plot original signal and continuous (without AR models)
plot_spectrum(continuous, 0);
ylim([-10 40]), title('Spectral analysis of continuous part')
plot_spectrum(z, 0);
ylim([-10 40]), title('Spectral analysis of original signal')