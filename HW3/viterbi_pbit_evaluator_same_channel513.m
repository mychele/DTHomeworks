%% This script performs simulation of BER for Viterbi
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec_viterbi = 5 : 15; % dB
L_data = 2.^[15 15 15 15 18 18 20 20 22 24 24] - 1;
snr_vec_viterbi = snr_vec_viterbi(1:9);
L_data = L_data(1:9);
if length(L_data) ~= length(snr_vec_viterbi), disp('Check L_data'), return, end

pbit_viterbi_same_ch = zeros(length(snr_vec_viterbi), 1);
n_biterr_viterbi_same_ch = zeros(length(snr_vec_viterbi), 1);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;

for snr_i = 1:length(snr_vec_viterbi)
    snr_ch = snr_vec_viterbi(snr_i);
    L_data_i = L_data(snr_i);
    
    % Load data
    load(strcat('inoutch', num2str(snr_ch), '.mat'));
    %packet = packet(1:25+(L_data_i-1)/2);
    %r = r(1:25+(L_data_i-1)/2+assumed_dly);
    % Sample to get r @ T
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    % Detection begins
    % Print progress update
    fprintf('Viterbi, snr = %d', snr_ch);
    % Compute
    [~, pbit_this, n_biterr_this] = viterbi(packet, x(1+assumed_dly-N1:end), hi, ...
        N1, N2, 0, N2, 25); % 25 is the length of the training sequence, that is only used
    % to train Viterbi and is not considered for pbit evaluation.
    pbit_viterbi_same_ch(snr_i) = pbit_this;
    n_biterr_viterbi_same_ch(snr_i) = n_biterr_this;
    fprintf('done!\n');
end

save('pbit_viterbi_same_ch513', 'pbit_viterbi_same_ch', 'n_biterr_viterbi_same_ch', 'snr_vec_viterbi');