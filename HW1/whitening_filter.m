%% Pass the signal through a whitening filter in order to equalize its spectrum
% The whitening filter will be obtained as the inverse of the AR model
% filter. This is needed to identify peaks. As we do not care about the
% phase in this stage, we can afford to pass it through something that is
% not with linear phase.

close all
clear all
clc

%% Import the signal
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length

%% Compute the AR model
[a, sigma_w] = arModel(3, autocorrelation(z, K/5));
[H, omega] = freqz(1, [1; a], K, 'whole');

%% Plot the two frequency responses
figure
plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2));
hold on
plot(omega/(2*pi), 10*log10((1/sigma_w)*abs(1./H).^2));
hold off
legend('AR filter', 'Inverse of AR filter');
title('Frequency response of the AR filter and of its inverse');
axis([0 1 -40 40]);

%% Plot the whitened result
figure
equalized_spectrum = abs((fft(z)).*(1./H)*(1/sigma_w));
plot((1:K-1)/K, 10*log10(equalized_spectrum(2:end))); hold on
title('Equalized spectrum of the given signal');
axis([0 1 -5 20]);