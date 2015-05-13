% This script solves the first problem.
clear
close all
clc
rng default % for reproducibility

Tc = 1;
T = 4*Tc;

%% Training sequence generation
L = 15;
Nseq = 10;
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
trainingseq = [mlseqdouble; mlseqdouble(1:2*Nseq)];
trainingsymbols = bitmap(trainingseq);

%% Generate the channel output

snr = 20; %dB
snr_lin = 10^(snr/10);
[r, sigma_w, q] = channel_output(trainingsymbols, T, Tc, snr_lin);


%% Choose N @T/4

maxN = 44;
error_func = zeros(maxN, 1);
for N_T4 = 1:maxN % N is the supposed length of the impulse response of the channel
    % Compute the supposed length of each branch
    n_short = mod(4-N_T4, 4); % Num branches with a shorter filter than others
    % N_i is the number of coefficients of the filter of the i-th branch.
    N_i(1:4-n_short) = ceil(N_T4/4);
    N_i(4-n_short + 1 : 4) = ceil(N_T4/4) - 1;
    x_for_ls = trainingsymbols(end-(L+max(N_i)-1) + 1 : end);
    d_for_ls = r(end - 4*(L+max(N_i)-1) + 1 - (length(q)-4): end - (length(q)-4));
    [~, r_hat] = h_estimation(x_for_ls, d_for_ls, L, N_i);
    d_no_trans = r(end-length(r_hat)+1- (length(q)-4) : end - (length(q)-4));
    error_func(N_T4) = sum(abs(r_hat - d_no_trans).^2)/length(r_hat); % TODO check if we should divide for L?
end

figure
plot(10*log10(error_func)), hold on, plot(10*log10(sigma_w*ones(1, maxN)))
title('Error functional estimating h @T/4')
xlabel('N_{T/4}'), ylabel('\epsilon [dB]')
grid on

% Fix N
N_T4 = 28;
N = ceil(N_T4 / 4);


%% Estimate q_hat @T/4 for the chosen N_T4

% Compute the supposed length of each branch
n_short = mod(4-N_T4, 4); % Num branches with a shorter filter than others
% N_i is the number of coefficients of the filter of the i-th branch.
N_i(1:4-n_short) = ceil(N_T4/4);
N_i(4-n_short + 1 : 4) = ceil(N_T4/4) - 1;
x_for_ls = trainingsymbols(end-(L+max(N_i)-1) + 1 : end);
d_for_ls = r(end - 4*(L+max(N_i)-1) + 1 - (length(q)-4): end - (length(q)-4));
[q_hat, ~] = h_estimation(x_for_ls, d_for_ls, L, N_i);

q_hat_vec = reshape(q_hat, numel(q_hat), 1);
figure, plot(0:length(q)-1, abs(q)), hold on, plot(0:length(q_hat_vec)-1, abs(q_hat_vec))
legend('q', 'qhat'), grid on, xlim([-0.2, length(q) + 0.2]), ylim([-0.05, max(abs(q_hat_vec) + 0.05)])
title('q, qhat for N = 28, estimated @T/4')
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
m_opt = m_opt - 1; % because of MATLAB indexing
init_offs = mod(m_opt, 4);



%% Choose N @T

maxN = 11;
error_func = zeros(maxN, 1);
for N = 1:maxN % N is the supposed length of the impulse response of the channel
    x_for_ls = trainingsymbols(end-(L+N-1) + 1 : end);
    d_for_ls = r(end - 4*(L+N-1) + 1 - (length(q)-4) + init_offs: 4 :end - (length(q)-4) + init_offs);
    [~, r_hat] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);    
    d_no_trans = d_for_ls(N : N+L-1);
    error_func(N) = sum(abs(r_hat - d_no_trans).^2)/length(r_hat);
end

figure
plot(10*log10(error_func)), hold on, plot(10*log10(sigma_w*ones(1, maxN)))
title('Error functional estimating h @T')
xlabel('N_{T}'), ylabel('\epsilon [dB]')
grid on

% Fix N
N = 7;


%% Estimate q_hat @T for the chosen N

x_for_ls = trainingsymbols(end-(L+N-1) + 1 : end);
d_for_ls = r(end - 4*(L+N-1) + 1 - (length(q)-4) + init_offs: 4 :end - (length(q)-4) + init_offs);
[q_hat, ~] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);

% Temp plot
figure, stem(init_offs:length(x_for_ls)-1+init_offs, real(x_for_ls))
hold on
plot(0:length(d_for_ls)-1, real(d_for_ls / q(11)))
title('d (normalized) and x, aligned according to t_0, @T')
legend('x', 'd'), grid on

q_hat_vec = q_hat.';
figure, plot(0:length(q)-1, abs(q)), hold on, stem(init_offs:4:4*length(q_hat_vec)-1+init_offs, abs(q_hat_vec))
legend('q', 'qhat'), grid on, xlim([-0.2, length(q) + 0.2]), ylim([-0.05, max(abs(q_hat_vec) + 0.05)])
title('q, qhat for N = 28, estimated @T')
xlabel('iT/4')


%% Interpolate q_hat to get h_i

% If estimation happened @T/4
% hi = q_hat_vec(init_offs+1:4:end);
% figure, plot(0:length(q)-1, abs(q)), hold on, plot(0:length(q_hat_vec)-1, abs(q_hat_vec)), hold on,
% stem(init_offs:4:length(q), abs(hi))
% legend('q', 'qhat', 'hi'), grid on, xlim([-0.2, length(q) + 0.2]), ylim([-0.05, max(abs(q_hat_vec) + 0.05)])
% title('q, qhat and hi')
% xlabel('iT/4')

% If estimation happened @T
hi = q_hat_vec;
figure, plot(0:length(q)-1, abs(q)), hold on, stem(init_offs:4:4*(length(hi)-1)+init_offs, abs(hi))
legend('q', 'qhat = hi'), grid on, xlim([-0.2, length(q) + 0.2]), ylim([-0.05, max(abs(q_hat_vec) + 0.05)])
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