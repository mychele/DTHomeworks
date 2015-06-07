%% OFDM BER
load('OFDM_uncoded_estimated.mat');
load('OFDM_coded_estimated.mat');
load('OFDM_coded_known.mat');
load('OFDM_uncoded_known.mat');

%% Plot
figure,
semilogy(snr_vec_coding_known, BER_coding_known, '-v'), hold on,
semilogy(snr_vec_nocoding_known, BER_nocoding_known, '-v'), hold on,
semilogy(snr_vec_coding_estimated, BER_coding_estimated, '--^'), hold on,
semilogy(snr_vec_nocoding_estimated, BER_nocoding_estimated, '--^')
legend('OFDM, known channel, coded', 'OFDM, known channel, uncoded', ...
    'OFDM, estimated channel, coded', 'OFDM, estimated channel, uncoded')
xlabel('SNR')
ylabel('Pbit')
ylim([10^-5, 10^-1])
xlim([0, 14])
grid on
title('BER for OFDM')

%% Comparison between channel estimations
load('OFDM_uncoded_estimated_2_long.mat');
load('OFDM_coded_estimated_2_long.mat');
load('OFDM_uncoded_estimated_long.mat');
load('OFDM_coded_estimated_long.mat');


%% Plot
figure,
semilogy(snr_vec_coding_estimated, mean(BER_coding_estimated, 2), '--.', 'MarkerSize', 3), hold on,
semilogy(snr_vec_nocoding_estimated, mean(BER_nocoding_estimated,2), '--.', 'MarkerSize', 3), hold on
semilogy(snr_vec_coding_estimated_2, mean(BER_coding_estimated_2, 2), '-v', 'MarkerSize', 3), hold on,
semilogy(snr_vec_nocoding_estimated_2, mean(BER_nocoding_estimated_2, 2), '-v', 'MarkerSize', 3)
legend('OFDM, estimated channel, coded, standard', 'OFDM, estimated channel, uncoded, standard', ...
    'OFDM, estimated channel, coded, new', 'OFDM, estimated channel, uncoded, new')
xlabel('SNR')
ylabel('Pbit')
ylim([10^-5, 10^-0.5])
grid on
title('BER for OFDM, comparison between the BER with the 2 different estimation methods')

%% Comparison between DFE and OFDM

load('Problem2_estch_coded.mat');
load('Problem2_knownch_coded.mat');
load('OFDM_coded_estimated.mat');
load('OFDM_coded_known.mat');

%% Plot
figure, 
semilogy(snr_vec_knownch_coded, Pbit_knownch_coded, '-v'), hold on
semilogy(snr_vec_estch_coded, Pbit_estch_coded, '--^'), hold on,
semilogy(snr_vec_coding_known, BER_coding_known, '-v'), hold on,
semilogy(snr_vec_coding_estimated, BER_coding_estimated, '--^'), hold on,
xlabel('\Gamma [dB]'), ylabel('Pbit'), grid on
ylim([10^-5, 10^-1])
xlim([0, 4])
legend('DFE, known channel', 'DFE, estimated channel', ...
    'OFDM, known channel', 'OFDM, estimated channel');
title('Comparison between Bit Error Rate for DFE and OFDM systems')
