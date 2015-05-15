% This script produces data that can be used by the receiver to perform
% detection. Use this to simulate tx with different SNRs.
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr = 14; % 6, 8, 10, 12, 14 % dB
L_data = 2^20-1;

%% Create, send and receive data with the given channel
[packet, r_T4, ~] = txrc(L_data, snr, T, Tc);

% estimate the channel using the first 100 samples (4*length(ts))
N = 7;
[ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);

% sample to get r @ T
init_offs = mod(m_opt, 4); % in T/4
t0 = floor(m_opt/4); % from now consider T = 1
T = 1;

%% Detection begins
rT = r_T4(init_offs+1:4:end); % data sampled in T
x = rT(t0+1:t0+1+length(packet)-1)/h(N1+1).'; % data normalized by h0, starting from t0
hi = h/h(N1+1).'; % impulse response normalized by h0

N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1 = N;   % FF filter: equal to the span of h
D = (M1-1);   % D is chosen large first and then decreased % (N-1)/2 + 2 or 3 for LE
M2 = N2 + M1 - 1 - D;      % FB filter: one less than the FF filter
[decisions, pbit, Jmin] = DFE_filter(packet, x, hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);


