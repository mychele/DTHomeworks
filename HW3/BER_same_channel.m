%% This script performs simulation of BER for LE, DFE, VA, FBA using
% the same estimate of the channel for the detection with each method
clear
close all
clc
rng default

snr_vec = 5 : 15; % dB
L_data = 2.^[15 15 15 15 18 18 20 20 22 23 24] - 1;
L_data_LE = 2.^[13 13 13 13 14 15 15 18 18 22 22] - 1;
L_data_DFE = 2.^[13 13 13 13 15 16 18 20 22 23 24] - 1;
if length(L_data) < length(snr_vec), disp('Check L_data'), return, end

pbit_LE_same_ch = zeros(length(snr_vec), 1);
n_biterr_LE_same_ch = zeros(length(snr_vec), 1);
pbit_DFE_same_ch = zeros(length(snr_vec), 1);
n_biterr_DFE_same_ch = zeros(length(snr_vec), 1);
pbit_viterbi_same_ch = zeros(length(snr_vec), 1);
n_biterr_viterbi_same_ch = zeros(length(snr_vec), 1);
pbit_fba_same_ch = zeros(length(snr_vec), 1);
n_biterr_fba_same_ch = zeros(length(snr_vec), 1);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;

for snr_i = 1:length(snr_vec)
    
    % Create, send and receive data, estimate channel and prepare for detection
    data_len = L_data(snr_i);
    snr_ch = snr_vec(snr_i);
    % Create, send and receive data with the given channel
    fprintf('Generating input symbols and channel output... ')
    [packet, r, sigma_w] = txrc(data_len, snr_ch, assumed_m_opt);
    fprintf('done!\n')
    
    % Estimate the channel using the first 25 samples (length(ts))
    fprintf('Estimating timing phase and IR... ')
    [ h, est_sigmaw ] = get_channel_info(r(assumed_dly+1 : 25+assumed_dly), N1, N2);
    fprintf('done!\n')
    
    % Normalize to get a QPSK constellation
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    % LE
    % Use only the first part of the data for LE and DFE, in order to speed
    % up the script (but enough to correctly estimate the BER)
    packet_LE = packet(1:25+(L_data_LE(snr_i)-1)/2);
    x_LE = x(1:25+(L_data_DFE(snr_i)-1)/2+assumed_dly);
    M1_le = 20;        % FF filter: equal to the span of h
    D_le = 15;
    M2_le = 0;      % FB filter not present in LE
    fprintf('LE for snr = %d ', snr_ch);
    [~, pbit, num_err, ~] = DFE_filter(packet_LE, x_LE(1+assumed_dly : assumed_dly+length(packet_LE)), ...
        hi, N1, N2, est_sigmaw, assumed_dly, D_le, M1_le, M2_le, 0);
    % DFE_filter with M2 = 0 acts as a LE
    pbit_LE_same_ch(snr_i, 1) = pbit;
    n_biterr_LE_same_ch(snr_i, 1) = num_err;
    fprintf('done!\n');
    
    % DFE
    packet_DFE = packet(1:25+(L_data_DFE(snr_i)-1)/2);
    x_DFE = x(1:25+(L_data_DFE(snr_i)-1)/2+assumed_dly);
    N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
    M1_dfe = 25;        % FF filter: equal to the span of h
    D_dfe = M1_dfe-1;
    M2_dfe = N2 + M1_dfe - 1 - D_dfe;  % FB filter
    fprintf('DFE for snr = %d ', snr_ch);
    [~, pbit, num_err, ~] = DFE_filter(packet_DFE, x_DFE(1+assumed_dly : assumed_dly+length(packet_DFE)), ...
        hi, N1, N2, est_sigmaw, assumed_dly, D_dfe, M1_dfe, M2_dfe, 0);
    pbit_LE_same_ch(snr_i, 1) = pbit;
    n_biterr_LE_same_ch(snr_i, 1) = num_err;
    fprintf('done!\n');
    
    % VA
    % Use the whole sequence
    fprintf('Viterbi, snr = %d ', snr_ch);
    [~, pbit, n_biterr] = viterbi(packet, x(1+assumed_dly-N1:end), hi, ...
        N1, N2, 0, N2, 25); % 25 is the length of the training sequence, that is only used
    % to train Viterbi and is not considered for pbit evaluation.
    pbit_viterbi_same_ch(snr_i) = pbit;
    n_biterr_viterbi_same_ch(snr_i) = n_biterr;
    fprintf('done!\n');
    
    % FBA
    if (snr_ch <= 13) % the BER for snr = 13 dB is already below 10^-5
        fprintf('FBA, snr = %d ', snr_ch);
        [~, pbit, n_biterr] = fba(packet, ...
            x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2);
        % 25 is the length of the training sequence, that is only used
        % to train Viterbi and is not considered for pbit evaluation.
        pbit_fba_same_ch(snr_i) = pbit;
        n_biterr_fba_same_ch(snr_i) = n_biterr;
        fprintf('done!\n');
    end
    
end

% Save the data computed, it will be processed by the BER_plot script
save('pbit_LE_same_channel', 'pbit_LE_same_ch', 'n_biterr_LE_same_ch', 'L_data_LE');
save('pbit_DFE_same_channel', 'pbit_DFE_same_ch', 'n_biterr_DFE_same_ch', 'L_data_DFE');
save('pbit_viterbi_same_channel', 'pbit_viterbi_same_ch', 'n_biterr_viterbi_same_ch', 'L_data');
save('pbit_fba_same_channel', 'pbit_fba_same_ch', 'n_biterr_fba_same_ch', 'L_data');

