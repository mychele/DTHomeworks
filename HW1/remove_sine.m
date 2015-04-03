% This script removes a sinusoidal interferer from a signal, using the lms
% algorithm.

% See page 224 for reference

%% Initialisation

% Clear stuff
close all;
clear all;
clc;

% Load actual signal we need to clean up
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length
autoc_z = autocorrelation(z, K/5);

N = 2;
upper_limit = 999; %MATLAB requires indices from 1 to 401
c = zeros(N, upper_limit + 1); % init c vector, no info -> set to 0
b = 5;
f0 = 0.771;
phi = pi;
n = 1:K;
x1 = b*cos(2*pi*f0*n + phi);
x2 = b*sin(2*pi*f0*n + phi);
e = zeros(1, upper_limit+1);
y = zeros(1, upper_limit);
d = z;

mu = 1/(autoc_z(1)*N); % actually mu must be > 0 and < 2/(N r_z(0))

for k = 1:upper_limit
    y(k) = x1(k) * c(1,k) + x2(k) * c(2,k);
    e(k) = d(k) - y(k);
    c(1, k + 1) = c(1, k) + mu*e(k)*conj(x1(k)); % update the filter, c(k+1) = c(k) + mu*e(k)*conj(z(k-1))
    c(2, k + 1) = c(2, k) + mu*e(k)*conj(x2(k)); % update the filter, c(k+1) = c(k) + mu*e(k)*conj(z(k-1))
end

% Plot c coefficients
for index = 1:N
    figure
    subplot(2, 1, 1)
    plot(1:upper_limit+1, real(c(index, :)))
    title(['Real part of c' int2str(index)]);
    subplot(2, 1, 2)
    plot(1:upper_limit+1, imag(c(index, :)))
    title(['Imaginary part of c' int2str(index)]);
end

% figure
% for iteration = 1:upper_limit
%     [H, omega] = freqz(1, c(:,iteration), K, 'whole');
%     plot(omega/2/pi, 10*log10(abs(H)))
%     pause(0.1)
% end

%figure, plot(1:upper_limit, 10*log10(abs(e).^2))
%title('Error function at each iteration');

figure,plot(10*log10(abs(fft(e))))