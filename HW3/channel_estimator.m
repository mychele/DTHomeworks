% This script solves the first problem.
clear
close all
clc
rng default % for reproducibility

Tc = 1;
T = 4*Tc;

%% Training sequence generation
L = 15;
maxN = 44;
maxNi = ceil(maxN/4);
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
trainingseq = [mlseqdouble; mlseqdouble(1:2*maxNi)];
trainingsymbols = bitmap(trainingseq);

%% Generate the channel output

snr = 20; %dB
snr_lin = 10^(snr/10);
[r, sigma_w, q] = channel_output(trainingsymbols, T, Tc, snr_lin);

%% Estimate qhat
error_func = zeros(maxN, 1);
for N = 1:maxN % N is the supposed length of the impulse response of the channel
    % Compute the supposed length of each branch
    n_short = mod(4-N, 4); % Num branches with a shorter filter than others
    % N_i is the number of coefficients of the filter of the i-th branch.
    N_i(1:4-n_short) = ceil(N/4);
    N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
    [~, r_hat] = h_estimation( trainingsymbols(end-(L+max(N_i)-1) + 1 : end), ...
        r(end - 4*(L+max(N_i)-1) + 1: end), L, N_i);
    d_no_trans = r(end-length(r_hat)+1 : end);
    error_func(N) = sum(abs(r_hat - d_no_trans).^2)/length(r_hat); % TODO check if we should divide for L?
end

figure
plot(10*log10(error_func)), hold on, plot(10*log10(sigma_w*ones(1, maxN)))
title('Error functional')
xlabel('N')
ylabel('\epsilon [dB]')

%% Plot q and qhat fot N = 28
N = 28;
% Compute the supposed length of each branch
n_short = mod(4-N, 4); % Num branches with a shorter filter than others
% N_i is the number of coefficients of the filter of the i-th branch.
N_i(1:4-n_short) = ceil(N/4);
N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
[q_hat, r_hat] = h_estimation( trainingsymbols(end-(L+max(N_i)-1) + 1 : end), ...
    r(end - 4*(L+max(N_i)-1) + 1: end), L, N_i);
d_no_trans = r(end-length(r_hat)+1 : end);
error_func(N) = sum(abs(r_hat - d_no_trans).^2)/length(r_hat); % TODO check if we should divide for L?

q_hat_vec = reshape(q_hat, numel(q_hat), 1);
figure, plot(0:length(q)-1, abs(q)), hold on, plot(0:length(q_hat_vec)-1, abs(q_hat_vec))
legend('q', 'qhat'), grid on, xlim([-0.2, length(q) + 0.2]), ylim([-0.05, max(abs(q_hat_vec) + 0.05)])
title('q, qhat for N = 28')
xlabel('iT/4')

%% Estimate t0 as in pag 617 bc
m_min = 0;
m_max = 43;
i = 1;
crossvec = zeros(m_max - m_min, 1);
for m = m_min:m_max
    r_part = r(m+1:T:m+1+4*(L-1)); % pick L samples from r, spaced by T
    crossvec(i) = abs(sum(r_part.*conj(trainingsymbols(1:L)))/L); % as in 7.269
    i = i + 1;
end

figure, plot(m_min:m_max, crossvec), grid on
ylabel('crosscorrelation, output sampled from m')
xlabel('m')
title('crosscorrelation between input & output')
xlim([-0.2, m_max + 0.2])

[~, m_opt] = max(crossvec);
m_opt = m_opt - 1; % beacuse of MATLAB indexing

%% Interpolate qhat to get h_i

init_offs = mod(m_opt, 4);
hi = q_hat_vec(init_offs+1:4:end);

figure, plot(0:length(q)-1, abs(q)), hold on, plot(0:length(q_hat_vec)-1, abs(q_hat_vec)), hold on,
stem(init_offs:4:length(q), abs(hi))
legend('q', 'qhat', 'hi'), grid on, xlim([-0.2, length(q) + 0.2]), ylim([-0.05, max(abs(q_hat_vec) + 0.05)])
title('q, qhat and hi')
xlabel('iT/4')

%% Possible way to define N1 and N2
[h0, i0] = max(hi); % in this way i0 correspond to t0 + 0*T
N1 = i0 - 1; % number of precursors
N2 = length(hi) - i0; % number of postcursors

%%%%% TODO: remember!!! Divide your received signal and everything for h0
%%%%% in order to get the same constellation of the tx at receiver

figure, stem(-N1:N2, abs(hi)), grid on, xlabel('i'), ylabel('hi'), title('hi vs i \in [-N1, N2]')
xlim([- N1 - 0.2, N2 + 0.2])