% This script produces data that can be used by the receiver to perform
% detection. Use this to simulate tx with different SNRs.
clear
close all
clc

Tc = 1;
T = 4 * Tc;
snr_vec = [6, 8, 10, 12, 14]; % dB
L_data = 2^20-1;

for i = 1:length(snr_vec)
    disp(i)
    snr = snr_vec(i);
    
    % create, send and receive data with the given channel
    [~, r_T4(:, i), ~] = txrc(L_data, snr, T, Tc);
    
    % estimate the channel using the first 100 samples (4*length(ts))
    N = 7;
    [ m_opt, h(:, i), est_sigmaw(i), N1, N2 ] = get_channel_info(r_T4(1:100, i), N, T);
    
    % sample to get r @ T
    % init_offs = mod(m_opt, 4); % in T/4
    % t0 = mod(m_opt,4); % from now consider T = 1
    % T = 1;
    % rT = r_T4(init_offs+1:4:end); % data sampled in T
    % x = rT/h(N1+1); % data normalized by h0
    % hi = h/h(N1+1); % impulse response normalized by h0
end
