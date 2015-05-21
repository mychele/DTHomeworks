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

pbit_viterbi = zeros(length(snr_vec_viterbi), numsim);
n_biterr_viterbi = zeros(length(snr_vec_viterbi), numsim);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;

% Useful data
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];
h = q(11:4:end);
h = h(:);
E_h = sum(abs(h).^2);
sigma_a_2 = 2;

numsim = 1;
%parpool(2);

for snr_i = 1:length(snr_vec_viterbi)
    thissnrstart = tic;
    L_data_curr = L_data(snr_i);
    snr_curr = snr_vec_viterbi(snr_i);
    for sim = 1:numsim     % perform sim_each simulations
        
        % --- Create, send and receive data, estimate channel and prepare for detection
        
        % Create, send and receive data with the given channel
        [packet, r, ~] = txrc(L_data_curr, snr_curr, assumed_m_opt);
        
        % Normalize x and h
        x = r / h(N1+1).';    % data normalized by h0
        hi = h / h(N1+1).';   % impulse response normalized by h0

        
        % --- Detection
        
        [~, pbit_this, n_biterr_this] = viterbi(packet, x(1+assumed_dly-N1:end), hi, ...
            N1, N2, 0, N2, 25); % 25 is the length of the training sequence, that is only used
                                % to train Viterbi and is not considered for pbit evaluation.
        pbit_viterbi(snr_i, sim) = pbit_this;
        n_biterr_viterbi(snr_i, sim) = n_biterr_this;
    end
    
    fprintf('SNR=%ddB completed in %.2f minutes\n', snr_curr, toc(thissnrstart)/60)
end

%delete(gcp);

save('pbit_viterbi_realch_1', 'pbit_viterbi', 'n_biterr_viterbi', 'snr_vec_viterbi');