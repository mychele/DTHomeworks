%% This script performs simulation of BER for Viterbi
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec_viterbi = 5 : 15; % dB
L_data = 2.^[18 20 20 20 20 20 20 20 22 22 22] - 1;
if length(L_data) ~= length(snr_vec_viterbi), disp('Check L_data'), return, end

sim_each = 4;
pbit_viterbi = zeros(length(snr_vec_viterbi), sim_each);
n_biterr_viterbi = zeros(length(snr_vec_viterbi), sim_each);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;


parpool(length(snr_vec_viterbi));

parfor snr_i = 1:length(snr_vec_viterbi)
    parforstart = tic;
    snr_curr = snr_vec_viterbi(snr_i);
    for sim = 1:sim_each     % perform sim_each simulations
        
        % --- Create, send and receive data, estimate channel and prepare for detection
        
        % Create, send and receive data with the given channel
        %fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data(snr_i), snr_curr, assumed_m_opt);
        %fprintf('done!\n')
        
        % Estimate the channel using the first 100 samples (4*length(ts))
        %fprintf('Estimating timing phase and IR... ')
        [ h, est_sigmaw ] = get_channel_info(r(assumed_dly+1:25+assumed_dly), N1, N2);
        %fprintf('done!\n')
        
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
    
    parforelapsed = toc(parforstart);
    fprintf('Completed %d simulations for SNR=%ddB in %.2f minutes\n', sim_each, snr_curr, parforelapsed/60)
end

delete(gcp);

save('pbit_viterbi', 'pbit_viterbi', 'n_biterr_viterbi', 'snr_vec_viterbi');