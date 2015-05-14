% This script produces data that can be used by the receiver to perform
% detection. Use this to simulate tx with different SNRs.
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec = [20];% 8, 10, 12, 14]; % dB
L_data = 127; %2^20-1;

for i = 1:length(snr_vec)
    disp(i)
    snr = snr_vec(i);
    
    % create, send and receive data with the given channel
    [packet, r_T4(:, i), ~] = txrc(L_data, snr, T, Tc);
    
    % estimate the channel using the first 100 samples (4*length(ts))
    N = 7;
    [ m_opt, h(:, i), est_sigmaw(i), N1, N2 ] = get_channel_info(r_T4(1:100, i), N, T);
    
    % sample to get r @ T
    init_offs = mod(m_opt, 4); % in T/4
    t0 = floor(m_opt/4); % from now consider T = 1
    %T = 1;
    
end
