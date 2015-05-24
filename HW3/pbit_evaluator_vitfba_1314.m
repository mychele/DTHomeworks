%% This script performs simulation of BER for Viterbi
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
L_data = 2.^[23] - 1;
snr_vec = 14;
if length(L_data) ~= length(snr_vec), disp('Check L_data'), return, end


% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;

% prepare data
for snr_i = 1:length(snr_vec)
    %% Create, send and receive data, estimate channel and prepare for detection
    
    data_len = L_data(snr_i);
    snr_ch = snr_vec(snr_i);
    % Create, send and receive data with the given channel
    fprintf('Generating input symbols and channel output... ')
    [packet, r, ~] = txrc(data_len, snr_ch, assumed_m_opt);
    fprintf('done!\n')
    
    % Estimate the channel using the first 100 samples (4*length(ts))
    fprintf('Estimating timing phase and IR... ')
    [ h, est_sigmaw ] = get_channel_info(r(assumed_dly+1 : 25+assumed_dly), N1, N2);
    fprintf('done!\n')
    
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    % Detection begins
    % Print progress update
    fprintf('Viterbi, snr = %d \n', snr_ch);
    % Compute
    [~, pbit_this_vit, n_biterr_this_vit] = viterbi(packet, x(1+assumed_dly-N1:end), hi, ...
        N1, N2, 0, N2, 25); % 25 is the length of the training sequence, that is only used
    % to train Viterbi and is not considered for pbit evaluation.
    fprintf('done!\n');
    
    save(strcat('pbit_viterbi_same_ch', num2str(snr_ch)), 'pbit_this_vit', 'n_biterr_this_vit', 'snr_ch');
    
    fprintf('FBA, snr = %d \n', snr_ch);
    [~, pbit_this_fba, n_biterr_this_fba] = fba(packet, ...
        x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2);
    fprintf('done!\n');
    
    save(strcat('pbit_fba_same_ch', num2str(snr_ch)), 'pbit_this_fba', 'n_biterr_this_fba', 'snr_ch');
end