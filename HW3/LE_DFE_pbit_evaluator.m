%% This script performs simulation of BER for LD, DFE
clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]; % dB

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
        %% Create, send and receive data with the given channel, just once for each snr
        T = 4*Tc;
        [packet, r_T4, ~] = txrc(L_data, snr_ch, T, Tc);
        
        % estimate the channel using the first 100 samples (4*length(ts))
        N = 7; % fixed, as specified by the teacher
        [ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);
        
        % sample to get r @ T
        init_offs = mod(m_opt, 4); % in T/4
        t0 = floor(m_opt/4); % from now consider T = 1
        T = 1;
        
        %% Detection begins
        rT = r_T4(init_offs+1:4:end); % data sampled in T
        x = rT(t0+1:t0+1+length(packet)-1)/h(N1+1).'; % data normalized by h0, starting from t0
        hi = h/h(N1+1).'; % impulse response normalized by h0
        
        N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
        M1 = 20;        % FF filter: equal to the span of h
        D = 14;
        M2 = 0;      % FB filter not present in LE
        % Print progress update
        fprintf('LE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
        % Compute only the functional
        [~, pbit, num_err, ~] = DFE_filter(packet, x, hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);
        pbitLE(snr_i, sim) = pbit;
        num_bit_errorLE(snr_i, sim) = num_err;
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
        %% Create, send and receive data with the given channel, just once for each snr
        T = 4*Tc;
        [packet, r_T4, ~] = txrc(L_data, snr_ch, T, Tc);
        
        % estimate the channel using the first 100 samples (4*length(ts))
        N = 7; % fixed, as specified by the teacher
        [ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);
        
        % sample to get r @ T
        init_offs = mod(m_opt, 4); % in T/4
        t0 = floor(m_opt/4); % from now consider T = 1
        T = 1;
        
        %% Detection begins
        rT = r_T4(init_offs+1:4:end); % data sampled in T
        x = rT(t0+1:t0+1+length(packet)-1)/h(N1+1).'; % data normalized by h0, starting from t0
        hi = h/h(N1+1).'; % impulse response normalized by h0
        
        N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
        M1 = 20;        % FF filter: equal to the span of h
        D = M1-1;
        M2 = N2 + M1 - 1 - D;  % FB filter
        % Print progress update
        % Print progress update
        fprintf('DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
        
        % Compute only the functional
        [~, pbit, num_err, ~] = DFE_filter(packet, x, hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);
        pbitDFE(snr_i, sim) = pbit;
        num_bit_errorDFE(snr_i, sim) = num_err;
    end
end

save('pbit_DFE', 'pbitDFE', 'num_bit_errorDFE');

delete(gcp);

%% Statistics

BER_LE = median(pbitLE, 2);
BER_DFE = median(pbitDFE, 2);

BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, BER_LE), hold on, semilogy(snr_vec, BER_DFE), hold on, semilogy(snr_vec, BER_ideal),
xlabel('snr [dB]'), ylabel('BER'), legend('LE, M1 = 20, D = 14', 'DFE, M1 = 20, D = 19', 'AWGN')
ylim([10^-6, 10^-1])

