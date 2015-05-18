%% This script performs simulation of BER for LE, DFE
clear
close all
clc
rng default

T = 1;
snr_vec = 6:2:14; %[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]; % dB
% From exercise 1
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;                   % t0 is @ T; TODO this is N1, we should refactor

%% LE pbit evaluation
L_data = 2^20-1; % be patient

sim_each = 64;
pbitLE = zeros(length(snr_vec), sim_each);
num_bit_errorLE = zeros(length(snr_vec), sim_each);

parpool(4);

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    parfor sim = 1:sim_each
        % perform sim_each simulations
        
        % Create, send and receive data with the given channel
        fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data, snr_ch, assumed_m_opt);
        fprintf('done!\n')
        
        % Estimate the channel using the first 100 samples (4*length(ts))
        fprintf('Estimating timing phase and IR... ')
        [ h, est_sigmaw, N1, N2 ] = get_channel_info(r(init_offs+1:25+init_offs), 0, 4, assumed_m_opt);
        fprintf('done!\n')
        
        % Sample to get r @ T
        x = r / h(N1+1).';         % data normalized by h0
        hi = h / h(N1+1).';         % impulse response normalized by h0
        
        % Detection begins
        N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
        M1 = 20;        % FF filter: equal to the span of h
        D = 16;
        M2 = 0;      % FB filter not present in LE
        % Print progress update
        fprintf('LE, snr = %d, M1 = %d, D = %d... ', snr_ch, M1, D);
        % Compute
        [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);
        pbitLE(snr_i, sim) = pbit;
        num_bit_errorLE(snr_i, sim) = num_err;
        fprintf('done!\n');
    end
end

save('pbit_LE', 'pbitLE', 'num_bit_errorLE');


%% DFE pbit evaluation

pbitDFE = zeros(length(snr_vec), sim_each);
num_bit_errorDFE = zeros(length(snr_vec), sim_each);

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    % perform sim_each simulations
    parfor sim = 1:sim_each
        % Create, send and receive data with the given channel
        fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data, snr_ch, assumed_m_opt);
        fprintf('done!\n')
        
        % Estimate the channel using the first 100 samples (4*length(ts))
        fprintf('Estimating timing phase and IR... ')
        [ h, est_sigmaw, N1, N2 ] = get_channel_info(r(init_offs+1:25+init_offs), 0, 4, assumed_m_opt);
        fprintf('done!\n')
        
        % Sample to get r @ T
        x = r / h(N1+1).';         % data normalized by h0
        hi = h / h(N1+1).';         % impulse response normalized by h0
        
        % Detection begins
        N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
        M1 = 25;        % FF filter: equal to the span of h
        D = M1-1;
        M2 = N2 + M1 - 1 - D;  % FB filter
        % Print progress update
        fprintf('DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
        % Compute
        [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);
        pbitDFE(snr_i, sim) = pbit;
        num_bit_errorDFE(snr_i, sim) = num_err;
        fprintf('done!\n');
    end
end

save('pbit_DFE', 'pbitDFE', 'num_bit_errorDFE');

delete(gcp);

%% Statistics

BER_LE = mean(pbitLE, 2);
BER_DFE = mean(pbitDFE, 2);

BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, BER_LE), hold on, semilogy(snr_vec, BER_DFE), hold on, semilogy(snr_vec, BER_ideal),
xlabel('snr [dB]'), ylabel('BER'), legend('LE, M1 = 20, D = 14', 'DFE, M1 = 20, D = 19', 'AWGN')
ylim([10^-6, 10^-1])

