% This script produces data that can be used by DFE
clear 
close all
clc

Tc = 1;
T = 4 * Tc;
snr = 20; % dB
L_data = 127;

% create, send and receive data with the given channel
[r_T4, ~] = txrc(L_data, snr, T, Tc);

% estimate the channel using the first 100 samples (4*length(ts))
N = 7;
[ m_opt, h_i, est_sigmaw, N1, N2] = get_channel_info(r_T4(1:100), N, T);

% sample to get r @ T
