%% This script performs simulation of BER for LE, DFE
clear
close all
clc
rng default
parpool(16);

T = 1;
snr_vec = 5 : 15; % dB
sim_each = 32;

% From exercise 1
assumed_dly = 2;
assumed_m_opt = 10;
N1 = 0;
N2 = 4;

%% LE pbit evaluation
L_data = 2.^[13 13 13 13 14 15 15 18 18 18 20] - 1;
if length(L_data) ~= length(snr_vec), disp('Check L_data'), return, end

pbitLE = zeros(length(snr_vec), sim_each);
num_bit_errorLE = zeros(length(snr_vec), sim_each);

% Detection begins
N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1_le = 20;        % FF filter: equal to the span of h
D_le = 15;
M2_le = 0;      % FB filter not present in LE

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    L_data_i = L_data(snr_i);
    parfor sim = 1:sim_each
        % Perform sim_each simulations
        
        % Create, send and receive data with the given channel
        fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data_i, snr_ch, assumed_m_opt);
        fprintf('done!\n')
        
        % Estimate the channel using the first 100 samples (4*length(ts))
        fprintf('Estimating timing phase and IR... ')
        [ h, est_sigmaw ] = get_channel_info(r(assumed_dly+1:25+assumed_dly), N1, N2);
        fprintf('done!\n')
        
        % Sample to get r @ T
        x = r / h(N1+1).';         % data normalized by h0
        hi = h / h(N1+1).';         % impulse response normalized by h0
        
        % Print progress update
        fprintf('LE, snr = %d, M1 = %d, D = %d... ', snr_ch, M1_le, D_le);
        % Compute
        [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
            hi, N1, N2, est_sigmaw, assumed_dly, D_le, M1_le, M2_le, 0);
        pbitLE(snr_i, sim) = pbit;
        num_bit_errorLE(snr_i, sim) = num_err;
        fprintf('done!\n');
    end
end

save('pbit_LE', 'pbitLE', 'num_bit_errorLE', 'snr_vec', 'M1_le', 'M2_le', 'D_le');

%% DFE pbit evaluation
L_data = 2.^[13 13 13 13 15 16 18 20 20 22 24] - 1;
if length(L_data) ~= length(snr_vec), disp('Check L_data'), return, end

pbitDFE = zeros(length(snr_vec), sim_each);
num_bit_errorDFE = zeros(length(snr_vec), sim_each);

% Detection begins
N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1_dfe = 25;        % FF filter: equal to the span of h
D_dfe = M1_dfe-1;
M2_dfe = N2 + M1_dfe - 1 - D_dfe;  % FB filter

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    L_data_i = L_data(snr_i);
    
    % perform sim_each simulations
    parfor sim = 1:sim_each
        % Create, send and receive data with the given channel
        fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data_i, snr_ch, assumed_m_opt);
        fprintf('done!\n')
        
        % Estimate the channel using the first 100 samples (4*length(ts))
        fprintf('Estimating timing phase and IR... ')
        [ h, est_sigmaw ] = get_channel_info(r(assumed_dly+1:25+assumed_dly), N1, N2);
        fprintf('done!\n')
        
        % Sample to get r @ T
        x = r / h(N1+1).';         % data normalized by h0
        hi = h / h(N1+1).';         % impulse response normalized by h0
        
        % Print progress update
        fprintf('DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1_dfe, D_dfe);
        % Compute
        [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
            hi, N1, N2, est_sigmaw, assumed_dly, D_dfe, M1_dfe, M2_dfe, 0);
        pbitDFE(snr_i, sim) = pbit;
        num_bit_errorDFE(snr_i, sim) = num_err;
        fprintf('done!\n');
    end
end

save('pbit_DFE', 'pbitDFE', 'num_bit_errorDFE', 'snr_vec', 'M1_dfe', 'M2_dfe', 'D_dfe');

delete(gcp);

%% Statistics

% BER_LE = mean(pbitLE, 2);
% BER_DFE = mean(pbitDFE, 2);
% 
% BER_ideal = BER_awgn(snr_vec);
% 
% figure, semilogy(snr_vec, BER_LE), hold on, semilogy(snr_vec, BER_DFE), hold on, semilogy(snr_vec, BER_ideal),
% xlabel('snr [dB]'), ylabel('BER'), legend('LE, M1 = 20, D = 14', 'DFE, M1 = 25, D = 24', 'AWGN')
% ylim([10^-6, 10^-1])

