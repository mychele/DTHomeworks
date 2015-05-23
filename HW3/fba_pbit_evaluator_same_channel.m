%% This script evaluates pbit for FBA using a fixed channel
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec_fba_same_channel = 5 : 13; % dB
L_data = 2.^[13 13 13 13 14 15 15 20 22] - 1;

if length(L_data) ~= length(snr_vec_fba_same_channel), disp('Check L_data'), return, end

pbit_fba_same_channel = zeros(length(snr_vec_fba_same_channel),1);
n_biterr_fba_same_channel = zeros(length(snr_vec_fba_same_channel),1);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;

for snr_i = 1:length(snr_vec_fba_same_channel)
    thissnrstart = tic;
    % --- Create, send and receive data, estimate channel and prepare for detection
    
    snr_ch = snr_vec_fba_same_channel(snr_i);
    L_data_i = L_data(snr_i);
    
    % Load data
    load(strcat('inoutch', num2str(snr_ch), '.mat'));
    packet = packet(1:25+(L_data_i-1)/2);
    r = r(1:25+(L_data_i-1)/2+assumed_dly);
    % Sample to get r @ T
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    % --- Detection
    [~, pbit_this, n_biterr_this] = fba(packet, ...
        x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2);
    % 25 is the length of the training sequence, that is only used
    % to train Viterbi and is not considered for pbit evaluation.
    pbit_fba_same_channel(snr_i) = pbit_this;
    n_biterr_fba_same_channel(snr_i) = n_biterr_this;
    
    
    fprintf('SNR=%ddB completed in %.2f minutes\n', snr_vec_fba_same_channel(snr_i), toc(thissnrstart)/60)
end

save('pbit_fba_same_channel', 'pbit_fba_same_channel', 'n_biterr_fba_same_channel', 'snr_vec_fba_same_channel');