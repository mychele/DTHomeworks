%% This script performs simulation of BER for LE, DFE
clear
close all
clc
rng default

T = 1;
snr_vec_1415 = 14 : 15; % dB
L_data = 2.^[24 25] - 1;
if length(L_data) ~= length(snr_vec_1415), disp('Check L_data'), return, end

% From exercise 1
assumed_dly = 2;
assumed_m_opt = 10;
N1 = 0;
N2 = 4;

%% LE pbit evaluation

pbitLE_same_channel = zeros(length(snr_vec_1415), 1);
num_bit_errorLE_same_channel = zeros(length(snr_vec_1415), 1);

pbitDFE_same_channel = zeros(length(snr_vec_1415), 1);
num_bit_errorDFE_same_channel = zeros(length(snr_vec_1415), 1);

for snr_i = 1:length(snr_vec_1415)
     %% Create, send and receive data, estimate channel and prepare for detection
    
    data_len = L_data(snr_i);
    snr_ch = snr_vec_1415(snr_i);
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
    
    % Detection begins
    N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
    M1_le = 20;        % FF filter: equal to the span of h
    D_le = 15;
    M2_le = 0;      % FB filter not present in LE
    % Print progress update
    fprintf('LE, snr = %d, M1 = %d, D = %d... ', snr_ch, M1_le, D_le);
    % Compute
    [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
        hi, N1, N2, est_sigmaw, assumed_dly, D_le, M1_le, M2_le, 0);
    pbitLE_same_channel(snr_i, 1) = pbit;
    num_bit_errorLE_same_channel(snr_i, 1) = num_err;
    fprintf('done!\n');
    
    save(strcat('pbit_LE_same_channel', num2str(snr_ch)), 'pbitLE_same_channel', 'num_bit_errorLE_same_channel', 'snr_ch', 'data_len');
    
    % Detection begins
    N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
    M1_dfe = 25;        % FF filter: equal to the span of h
    D_dfe = M1_dfe-1;
    M2_dfe = N2 + M1_dfe - 1 - D_dfe;  % FB filter
    % Print progress update
    fprintf('DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1_dfe, D_dfe);
    % Compute
    [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
        hi, N1, N2, est_sigmaw, assumed_dly, D_dfe, M1_dfe, M2_dfe, 0);
    pbitDFE_same_channel(snr_i, 1) = pbit;
    num_bit_errorDFE_same_channel(snr_i, 1) = num_err;
    fprintf('done!\n');
    
    save(strcat('pbit_DFE_same_channel', num2str(snr_ch)), 'pbitDFE_same_channel', 'num_bit_errorDFE_same_channel', 'snr_ch', 'data_len');


end

% %% DFE pbit evaluation
% L_data = 2.^[12 13 13 13 15 16 18 20 22] - 1;
% 
% pbitDFE_same_channel = zeros(length(snr_vec_1415), 1);
% num_bit_errorDFE_same_channel = zeros(length(snr_vec_1415), 1);
% 
% for snr_i = 1:length(snr_vec_1415)
%     snr_ch = snr_vec_1415(snr_i);
%     L_data_i = L_data(snr_i);
%     
%     % Load data
%     load(strcat('inoutch', num2str(snr_ch), '.mat'));
%     packet = packet(1:25+(L_data_i-1)/2);
%     r = r(1:25+(L_data_i-1)/2+assumed_dly);
%     
%     % Sample to get r @ T
%     x = r / h(N1+1).';         % data normalized by h0
%     hi = h / h(N1+1).';         % impulse response normalized by h0
%     
%     % Detection begins
%     N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
%     M1_dfe = 25;        % FF filter: equal to the span of h
%     D_dfe = M1_dfe-1;
%     M2_dfe = N2 + M1_dfe - 1 - D_dfe;  % FB filter
%     % Print progress update
%     fprintf('DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1_dfe, D_dfe);
%     % Compute
%     [~, pbit, num_err, ~] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), ...
%         hi, N1, N2, est_sigmaw, assumed_dly, D_dfe, M1_dfe, M2_dfe, 0);
%     pbitDFE_same_channel(snr_i, 1) = pbit;
%     num_bit_errorDFE_same_channel(snr_i, 1) = num_err;
%     fprintf('done!\n');
% end
% 
% save('pbit_DFE_same_channel', 'pbitDFE_same_channel', 'num_bit_errorDFE_same_channel', 'snr_vec_513');
% 
% %% Statistics
% 
% % BER_LE = mean(pbitLE_same_channel, 2);
% % BER_DFE = mean(pbitDFE_real_channel, 2);
% %
% BER_ideal = BER_awgn(snr_vec_1415);
% %
% figure, semilogy(snr_vec_1415, pbitLE_same_channel), hold on, semilogy(snr_vec_1415, pbitDFE_same_channel), hold on, semilogy(snr_vec_1415, BER_ideal),
% xlabel('snr [dB]'), ylabel('BER'), legend('LE, M1 = 20, D = 14', 'DFE, M1 = 20, D = 19', 'AWGN')
% ylim([10^-6, 10^-1])

