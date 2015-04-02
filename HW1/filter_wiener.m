close all;
%clear all;
clc;

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length
autoc = autocorrelation(z, length(z)/5);

%% High pass filter

z = filter(hp30, 1, z);

plot_spectrum(z);

%% Filter pg 137 BC
N = 30;
w0 = 2*pi*0.78; % 0.78 is probably our freq
vec =  w0 * (0:N-1);
E_w0 = exp(1i * vec).';

A = sqrt(30); % see plot of real autoc, rule of thumb
B=A;
lambda = A^2/autoc(1); % autoc(1) should be sigma^2 of the noise, if white

gain = B/A;

D = 100; % delay

c_opt = gain * lambda * exp(- 1i * w0 * D) * E_w0 / (1 + N*lambda);

[H, w] = freqz(c_opt, 1, 'whole');

figure, plot(w/(2*pi), 20*log10(H))
axis([0, 1, -40, 5])
figure, plot(w/(2*pi), angle(H))

% Filter using the Weiner filter
spectra = filter(c_opt, 1, z);

% Plot original and filtered signal
figure
subplot(2,2,1)
plot(real(z))
title('Real part of original signal');
subplot(2,2,3)
plot(imag(z), 'r')
title('Imaginary part of original signal');
subplot(2,2,2)
plot(real(spectra))
title('Real part of spectral line');
subplot(2,2,4)
plot(imag(spectra), 'r')
title('Imaginary part of spectral line');

plot_spectrum(spectra);

N_corr = length(spectra)/5;
autoc = autocorrelation(spectra, N_corr);
% useful only to know which is the knee of the sigma_w
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
    [~, sigma_w(N)] = arModel(N, autoc);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))% Uncomment to plot sigma

%the knee is apparently at N = 2
%compute the vector of coefficients a
N = 2;
[a, sigma_w] = arModel(N, autoc);
[H, omega] = freqz(1, [1; a], K, 'whole');

figure
plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
axis([0, 1, -40, 10])


