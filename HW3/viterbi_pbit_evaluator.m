%% This script performs simulation of BER for LD, DFE
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec = 6 : 2 : 14; % dB
L_data = 2.^ [18 20 20 20 20] - 1;
if length(L_data) ~= length(snr_vec), disp('Check L_data'), return, end

sim_each = 8;
pbit_viterbi = zeros(length(snr_vec), sim_each);
n_biterr_viterbi = zeros(length(snr_vec), sim_each);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;


parpool(8);

for snr_i = 1:length(snr_vec)
    snr_curr = snr_vec(snr_i);
    parfor sim = 1:sim_each     % perform sim_each simulations
        
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
        
        [~, pbit_this, n_biterr_this] = viterbi(packet, [x(1+assumed_dly-N1:end); zeros(assumed_dly-N1, 1)], hi, N1, N2, 0, N2, 25);
        pbit_viterbi(snr_i, sim) = pbit_this;
        n_biterr_viterbi(snr_i, sim) = n_biterr_this;
    end
    
    fprintf('%.1f%% completed\n', snr_i*100/length(snr_vec))
end

delete(gcp);

save('pbit_viterbi', 'pbit_viterbi', 'n_biterr_viterbi');

%% Statistics

BER_viterbi = median(pbit_viterbi, 2);
BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, BER_viterbi), hold on, semilogy(snr_vec, BER_ideal)
xlabel('snr [dB]'), ylabel('BER'), legend('Viterbi', 'AWGN')
ylim([10^-6, 10^-1])