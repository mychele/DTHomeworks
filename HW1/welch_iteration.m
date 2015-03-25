close all;
clear all;
clc;

% pg 87 Benvenuto Cherubini - iterative improvement over Welch estimate

% load data
z = load('data for hw1.mat');
% make a column vector
z1 = z.z.';

z1 = z1 - mean(z1);

K = length(z1);

M = [];
for D = 100:100:500      % window size

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
    M = [M, sum(P_per_w, 2)/N_s];
end

plot((1:K)/K, 10*log10(M))
xlabel('Normalized frequency');
ylabel('PSD estimate [dB]');
axis([0 1 -5 40])


