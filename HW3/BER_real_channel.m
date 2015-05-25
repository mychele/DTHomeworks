%% This script performs simulation of BER for LE, DFE, VA, FBA using
% the real channel for the detection with each method
clear
close all
clc
rng default

% Given channel
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];
h = q(11:4:end);
h = h(:);
E_h = sum(abs(h).^2);
sigma_a_2 = 2;

snr_vec = 5 : 15; % dB
L_data = 2.^[15 15 15 15 18 18 20 20 22 23 23] - 1;
if length(L_data) < length(snr_vec), disp('Check L_data'), return, end

pbit_LE_real_ch = zeros(length(snr_vec), 1);
n_biterr_LE_real_ch = zeros(length(snr_vec), 1);
pbit_DFE_real_ch = zeros(length(snr_vec), 1);
n_biterr_DFE_real_ch = zeros(length(snr_vec), 1);
pbit_viterbi_real_ch = zeros(length(snr_vec), 1);
n_biterr_viterbi_real_ch = zeros(length(snr_vec), 1);
pbit_fba_real_ch = zeros(length(snr_vec), 1);
n_biterr_fba_real_ch = zeros(length(snr_vec), 1);

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
    
    % Normalize to get a QPSK constellation
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    % LE
    M1_le = 20;        % FF filter: equal to the span of h
    D_le = 15;
    M2_le = 0;      % FB filter not present in LE
    fprintf('LE for snr = %d ', snr_ch);
    [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
        hi, N1, N2, sigma_w, assumed_dly, D_le, M1_le, M2_le, 0);
    % DFE_filter with M2 = 0 acts as a LE
    pbit_LE_real_ch(snr_i, 1) = pbit;
    n_biterr_LE_real_ch(snr_i, 1) = num_err;
    fprintf('done!\n');
    
    % DFE
    N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
    M1_dfe = 25;        % FF filter: equal to the span of h
    D_dfe = M1_dfe-1;
    M2_dfe = N2 + M1_dfe - 1 - D_dfe;  % FB filter
    fprintf('DFE for snr = %d ', snr_ch);
    [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
        hi, N1, N2, sigma_w, assumed_dly, D_dfe, M1_dfe, M2_dfe, 0);
    pbit_LE_real_ch(snr_i, 1) = pbit;
    n_biterr_LE_real_ch(snr_i, 1) = num_err;
    fprintf('done!\n');
    
    if(snr_ch <= 13) % the BER for snr = 13 dB is already below 10^-5
        % VA
        % Use the whole sequence
        fprintf('Viterbi, snr = %d ', snr_ch);
        [~, pbit, n_biterr] = viterbi(packet, x(1+assumed_dly-N1:end), hi, ...
            N1, N2, 0, N2, 25); % 25 is the length of the training sequence, that is only used
        % to train Viterbi and is not considered for pbit evaluation.
        pbit_viterbi_real_ch(snr_i) = pbit;
        n_biterr_viterbi_real_ch(snr_i) = n_biterr;
        fprintf('done!\n');
        
        % FBA
        fprintf('FBA, snr = %d ', snr_ch);
        [~, pbit, n_biterr] = fba(packet, ...
            x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2);
        % 25 is the length of the training sequence, that is only used
        % to train Viterbi and is not considered for pbit evaluation.
        pbit_fba_real_ch(snr_i) = pbit;
        n_biterr_fba_real_ch(snr_i) = n_biterr;
        fprintf('done!\n');
    end
end

% Save the data computed, it will be processed by the BER_plot script
save('pbit_LE_real_channel', 'pbit_LE_real_ch', 'n_biterr_LE_real_ch', 'L_data_LE');
save('pbit_DFE_real_channel', 'pbit_DFE_real_ch', 'n_biterr_DFE_real_ch', 'L_data_DFE');
save('pbit_viterbi_real_channel', 'pbit_viterbi_real_ch', 'n_biterr_viterbi_real_ch', 'L_data');
save('pbit_fba_real_channel', 'pbit_fba_real_ch', 'n_biterr_fba_real_ch', 'L_data');

