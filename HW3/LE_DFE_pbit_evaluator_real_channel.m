%% This script performs simulation of BER for LE, DFE
clear
close all
clc
rng default
parpool(8);

T = 1;
snr_vec = 5 : 15; % dB
L_data = 2.^[18 20 20 20 20 20 20 20 22 22 22] - 1;
if length(L_data) ~= length(snr_vec), disp('Check L_data'), return, end
sim_each = 1;

% From exercise 1
assumed_dly = 2;
assumed_m_opt = 10;
N1 = 0;
N2 = 4;

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


%% LE pbit evaluation

pbitLE_real_channel = zeros(length(snr_vec), sim_each);
num_bit_errorLE_real_channel = zeros(length(snr_vec), sim_each);

parfor snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    L_data_i = L_data(snr_i);
    for sim = 1:sim_each
        % Perform sim_each simulations
        
        % Create, send and receive data with the given channel
        fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data_i, snr_ch, assumed_m_opt);
        fprintf('done!\n')
        
        % Sample to get r @ T
        x = r / h(N1+1).';         % data normalized by h0
        hi = h / h(N1+1).';         % impulse response normalized by h0
        
        % Detection begins
        N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
        M1 = 20;        % FF filter: equal to the span of h
        D = 15;
        M2 = 0;      % FB filter not present in LE
        % Print progress update
        fprintf('LE, snr = %d, M1 = %d, D = %d... ', snr_ch, M1, D);
        % Compute
        [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
                                    hi, N1, N2, sigma_w, assumed_dly, D, M1, M2, 0);
        pbitLE_real_channel(snr_i, sim) = pbit;
        num_bit_errorLE_real_channel(snr_i, sim) = num_err;
        fprintf('done!\n');
    end
end

save('pbit_LE_real_channel', 'pbitLE_real_channel', 'num_bit_errorLE_real_channel', 'snr_vec');

%% DFE pbit evaluation

pbitDFE_real_channel = zeros(length(snr_vec), sim_each);
num_bit_errorDFE_real_channel = zeros(length(snr_vec), sim_each);

parfor snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    L_data_i = L_data(snr_i);

    % perform sim_each simulations
    for sim = 1:sim_each
        % Create, send and receive data with the given channel
        fprintf('Generating input symbols and channel output... ')
        [packet, r, sigma_w] = txrc(L_data_i, snr_ch, assumed_m_opt);
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
        [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ... 
                                        hi, N1, N2, sigma_w, assumed_dly, D, M1, M2, 0);
        pbitDFE_real_channel(snr_i, sim) = pbit;
        num_bit_errorDFE_real_channel(snr_i, sim) = num_err;
        fprintf('done!\n');
    end
end

save('pbit_DFE_real_channel', 'pbitDFE_real_channel', 'num_bit_errorDFE_real_channel', 'snr_vec');

delete(gcp);

%% Statistics

BER_LE = mean(pbitLE_real_channel, 2);
BER_DFE = mean(pbitDFE_real_channel, 2);

BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, BER_LE), hold on, semilogy(snr_vec, BER_DFE), hold on, semilogy(snr_vec, BER_ideal),
xlabel('snr [dB]'), ylabel('BER'), legend('LE, M1 = 20, D = 14', 'DFE, M1 = 20, D = 19', 'AWGN')
ylim([10^-6, 10^-1])

