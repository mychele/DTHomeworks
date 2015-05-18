%% This script performs simulation of BER for LD, DFE
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec = 6 : 2 : 14; % dB
L_data = 2.^ [15, 15, 18, 20, 20] - 1;
if length(L_data) ~= length(snr_vec), disp('Check L_data'), return, end

sim_each = 10;
pbit = zeros(length(snr_vec), sim_each);
n_biterr = zeros(length(snr_vec), sim_each);

%parpool(2);

for snr_i = 1:length(snr_vec)
    snr_curr = snr_vec(snr_i);
    for sim = 1:sim_each     % perform sim_each simulations
        
        % --- Create, send and receive data, estimate channel and prepare for detection
        
        % Create, send and receive data with the given channel
        %fprintf('Generating input symbols and channel output... ')
        [packet, r_T4, ~] = txrc(L_data(snr_i), snr_curr, T, Tc);
        %fprintf('done!\n')
        
        % Estimate the channel using the first 100 samples (4*length(ts))
        N1 = 0;
        N2 = 4;
        assumed_dly = 2;
        m_opt = 10; % these 4 quantities are assumed from ex 1
        %fprintf('Estimating timing phase and IR... ')
        [ h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100+T*assumed_dly), N1, N2, m_opt, T);
        %fprintf('done!\n')
        
        % Sample to get r @ T
        init_offs = mod(m_opt, T);  % offset in T/4
        t0 = N1;                    % t0 is @ T; TODO this is N1, we should refactor
        rT = r_T4(init_offs+1:T:end); % data sampled in T
        x = rT / h(N1+1).';         % data normalized by h0
        hi = h / h(N1+1).';         % impulse response normalized by h0
        
        % --- Detection
        
        [~, pbit_this, n_biterr_this] = viterbi(packet(1:end-assumed_dly), x(1+assumed_dly:end), hi, N1, N2, 0, N2);
        pbit(snr_i, sim) = pbit_this;
        n_biterr(snr_i, sim) = n_biterr_this;
    end
    
    fprintf('%.1f%% completed\n', snr_i*100/length(snr_vec))
end

%delete(gcp);

save('BER_viterbi', 'pbit', 'n_biterr');

%% Statistics

BER_viterbi = median(pbit, 2);

BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, BER_viterbi), hold on, semilogy(snr_vec, BER_ideal)
xlabel('snr [dB]'), ylabel('BER'), legend('Viterbi', 'AWGN')
ylim([10^-6, 10^-1])