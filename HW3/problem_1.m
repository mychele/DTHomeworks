% This script solves the first problem.
clear, clc, close all

%rng default % for reproducibility

Tc = 1;
T = 4*Tc;
sigma_a = 2; % for the QPSK

%% Training sequence generation

L = 15;
Nseq = 10;
trainingsymbols = ts_generation(L, Nseq);

%% Generate the channel output

snr = 20; %dB
snr_lin = 10^(snr/10);
[r, sigma_w, q] = channel_output(trainingsymbols, T, Tc, snr_lin);

%% Estimate t0 as in pag 617 bc

m_min = 0;
m_max = 40;
crossvec = zeros(m_max - m_min, 1);
for m = m_min : m_max
    r_part = r(m+1 : T : m+1+T*(L-1));  % pick L samples from r, spaced by T
    crossvec(m+1) = abs(sum(r_part.*conj(trainingsymbols(1:L)))/L); % as in 7.269
end

[~, m_opt] = max(crossvec);
m_opt = m_opt - 1; % because of MATLAB indexing
init_offs = mod(m_opt, T);
delay = floor(m_opt / T);   % timing phase @T, that is the delay of the channel @T

%% Error functional for different N1 and N2 (@T), given the delay

maxN1 = delay; % we already know it can't be larger than this
maxN = 11;
error_func = zeros(maxN, maxN1);
rT = r(init_offs+1 : T : end);
for N1 = 0:maxN1
    for N2 = 0:maxN-N1-1
        N = N1+N2+1;
        x_for_ls = trainingsymbols(1 : L+N-1);
        d_for_ls = rT(1+delay-N1 : (delay-N1)+(L+N-1));
        %d_for_ls = r(end - T*(L+N-1) + 1 - (length(q)-4) + init_offs + T*(delay-N1): T :end - (length(q)-4) + T*(delay-N1) + init_offs);
        % Now d_for_ls is delayed by N1 samples wrt x_for_ls.
        [~, r_hat] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);
        % We discarded the first t0-N1 samples, so the "perceived delay" is N1 < t0.
        % We estimated the channel with N coefficients disregarding the first t0-N1,
        % and the IR we get is a version of the actual one shifted left by t0-N1.
        d_no_trans = d_for_ls(N : N+L-1);
        error_func(N-N1, N1+1) = sum(abs(r_hat - d_no_trans).^2)/length(r_hat);
    end
end

figure
plot(0:maxN-1, 10*log10(error_func))
title('Error functional estimating h @T, given t_0 and varying N_1')
legend('N_1 = 0', 'N_1 = 1', 'N_1 = 2')
xlabel('N_2'), ylabel('\epsilon [dB]'), grid on, xlim([0 8])


%% Channel estimation, and computation of Lambda_n, for two choices of N1 and N2

% --- N1 = 0

N1 = 0; N2 = 4;
N = N1+N2+1;
x_for_ls = trainingsymbols(1 : L+N-1);
d_for_ls = rT(1+delay-N1 : (delay-N1)+(L+N-1));
% Now d_for_ls is delayed by N1 samples wrt x_for_ls.
[hi_0, r_hat] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);
d_no_trans = d_for_ls(N : N+L-1);
est_sigmaw = sum(abs(r_hat - d_no_trans).^2)/length(r_hat);

h_true = q(1+init_offs + T*(delay-N1) : T : end);
lambdan_0_estim = est_sigmaw / (sigma_a * sum(abs(hi_0 - h_true).^2));
lambdan_0_true = (L+1) * (L+1-N) / (N * (L+2-N));


% --- N1 = 2

N1 = 2; N2 = 4;
N = N1+N2+1;
x_for_ls = trainingsymbols(1 : L+N-1);
d_for_ls = rT(1+delay-N1 : (delay-N1)+(L+N-1));
% Now d_for_ls is delayed by N1 samples wrt x_for_ls.
[hi_2, r_hat] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);
d_no_trans = d_for_ls(N : N+L-1);
est_sigmaw = sum(abs(r_hat - d_no_trans).^2)/length(r_hat);

h_true = q(1+init_offs + T*(delay-N1) : T : end);
lambdan_2_estim = est_sigmaw / (sigma_a * sum(abs(hi_2 - h_true).^2));
lambdan_2_true = (L+1) * (L+1-N) / (N * (L+2-N));


% --- Output Lambda_n

fprintf('N1 = 2\n')
fprintf('Lambda_n estimate: %.2f\n', lambdan_2_estim)
fprintf('Lambda_n true:     %.2f\n', lambdan_2_true)
fprintf('N1 = 0\n')
fprintf('Lambda_n estimate: %.2f\n', lambdan_0_estim)
fprintf('Lambda_n true:     %.2f\n', lambdan_0_true)



%% Plot hi for two choices of N1 and N2

figure, hold on
stem(-2:N2, abs(hi_2), 'x', 'Color', [0.85 0.325 0.098])
stem(0:N2, abs(hi_0), 'o')
grid on, box on
xlabel('i'), ylabel('|h_i|'), title('Estimated |h_i| vs i')
xlim([-N1 - 0.2, N2 + 0.2]), legend('|h_i| with N_1=0', '|h_i| with N_1=2')