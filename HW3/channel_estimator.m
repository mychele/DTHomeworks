% This script solves the first problem.
clear
close all
clc
rng default % for reproducibility

Tc = 1;
T = 4*Tc;

%% Training sequence generation
L = 15;
N = 10;
mlseq = MLsequence(L);

% Replace every 0 with two 0s and every 1 with two 1s to put it into the
% bitmap
mlseqdouble = zeros(2*L,1);
for i = 1:L
    switch mlseq(i)
        case 0
            mlseqdouble(2*i-1) = 0;
            mlseqdouble(2*i) = 0;
        case 1
            mlseqdouble(2*i-1) = 1;
            mlseqdouble(2*i) = 1;
    end
end

% Repeat the sequence and bitmap it to get the symbols
trainingseq = [mlseqdouble; mlseqdouble(1:2*N)];
trainingsymbols = bitmap(trainingseq);

%% Generate the channel output

snr = 20; %dB
snr_lin = 10^(snr/10);
[r, sigma_w] = channel_output(trainingsymbols, T, Tc, snr_lin);

%% Estimate qhat
maxN = 44;
error_func = zeros(maxN, 1);
for N = 1:maxN % N is the supposed length of the impulse response of the channel
    % Compute the supposed length of each branch
    n_short = mod(4-N, 4); % Num branches with a shorter filter than others
    % N_i is the number of coefficients of the filter of the i-th branch.
    N_i(1:4-n_short) = ceil(N/4);
    N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
    [h_hat, d_hat] = h_estimation( trainingsymbols(end-(L+max(N_i)-1) + 1 : end), ...
        r(end - 4*(L+max(N_i)-1) + 1: end), L, N_i);
    d_no_trans = r(end-length(d_hat)+1 : end);
    error_func(N) = sum(abs(d_hat - d_no_trans).^2)/length(d_hat); % TODO check if we should divide for L?
end

figure
plot(10*log10(error_func)), hold on, plot(10*log10(sigma_w*ones(1, maxN)))
