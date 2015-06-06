%% Load OFDM BER
load('OFDM_uncoded_estimated.mat');
load('OFDM_coded_estimated.mat');
load('OFDM_coded_known.mat');
load('OFDM_uncoded_known.mat');

%% Plot
figure,
semilogy(snr_vec_coding_known, BER_coding_known), hold on,
semilogy(snr_vec_nocoding_known, BER_nocoding_known), hold on,
semilogy(snr_vec_coding_estimated, BER_coding_estimated(:, 2), '--'), hold on,
semilogy(snr_vec_nocoding_estimated, BER_nocoding_estimated(:, 7), '--')
legend('OFDM, known channel, coded', 'OFDM, known channel, uncoded', ...
    'OFDM, estimated channel, coded', 'OFDM, estimated channel, uncoded')
xlabel('SNR')
ylabel('BER')
ylim([10^-5, 10^-0.5])
grid on
title('BER for OFDM')

%% New estimation method
load('OFDM_uncoded_estimated_2.mat');
load('OFDM_coded_estimated_2.mat');

%% Plot
figure,
semilogy(snr_vec_coding_estimated, median(BER_coding_estimated, 2), '--'), hold on,
semilogy(snr_vec_nocoding_estimated, BER_nocoding_estimated(:, 7), '--'), hold on
semilogy(snr_vec_coding_estimated_2, median(BER_coding_estimated_2, 2)), hold on,
semilogy(snr_vec_nocoding_estimated_2, BER_nocoding_estimated_2(:, 7), '--')
legend('OFDM, estimated channel, coded, standard', 'OFDM, estimated channel, uncoded, standard', ...
    'OFDM, estimated channel, coded', 'OFDM, estimated channel, uncoded')
xlabel('SNR')
ylabel('BER')
ylim([10^-5, 10^-0.5])
grid on
title('BER for OFDM')
