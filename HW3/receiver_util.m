% This script produces data that can be used by DFE
clear 
close all
clc

Tc = 1;
T = 4 * Tc;
snr = 20; % dB
L_data = 127;

% create, send and receive data with the given channel
[packet, r_T4, ~] = txrc(L_data, snr, T, Tc);

% estimate the channel using the first 100 samples (4*length(ts))
N = 7;
[ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);

% sample to get r @ T
init_offs = mod(m_opt, 4); % in T/4
t0 = mod(m_opt,4); % from now consider T = 1
T = 1;
rT = r_T4(init_offs+1:4:end); % data sampled in T
x = rT/h(N1+1); % data normalized by h0
hi = h/h(N1+1); % impulse response normalized by h0
