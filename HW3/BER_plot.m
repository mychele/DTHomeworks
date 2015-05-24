clear, clc, close all

%% Group all the BER curves

load('pbit_DFE.mat');
load('pbit_LE.mat');
load('pbit_viterbi.mat');
load('BER_awgn_sim.mat');
load('pbit_fba_513.mat');
load('pbit_fba_1415.mat');
pbit_fba = [pbit_fba_513(:, 2); pbit_fba_1415];

%% Statistics for estimated channel, with only one realization
snr_vec = 5:15;

figure
semilogy(snr_vec, pbitLE(:, 3)), hold on,
semilogy(snr_vec, pbitDFE(:, 2))
semilogy(snr_vec_viterbi, pbit_viterbi(:, 3))
semilogy(snr_vec, pbit_fba)
semilogy(snr_vec, BER_awgn(snr_vec))
semilogy(snr_vec, pbit_AWGN_sim)
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'Viterbi', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-5.5, 10^-1]), grid on
title('BER for one realization, estimated channel');


%% Statistics for real channel

load('pbit_LE_real_channel.mat');
load('pbit_DFE_real_channel.mat');
load('pbit_fba_real_channel.mat');
load('pbit_viterbi_realch.mat')

figure
semilogy(snr_vec, pbitLE_real_channel), hold on
semilogy(snr_vec, pbitDFE_real_channel)
semilogy(snr_vec_viterbi, pbit_viterbi)
semilogy(snr_vec_fba_real_channel, pbit_fba_real_channel)
semilogy(snr_vec, BER_awgn(snr_vec))
semilogy(snr_vec, pbit_AWGN_sim)
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'Viterbi', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-7, 10^-1]), xlim([5, 15]), grid on
title('BER for methods using the actual channel impulse response');

%% Statistics for same estimate

snr_vec = 5:13; 

load('pbit_LE_same_channel.mat');
load('pbit_DFE_same_channel.mat');
load('pbit_viterbi_same_channel.mat');
load('pbit_fba_same_channel.mat');

figure
semilogy(snr_vec, pbitLE_same_channel), hold on,
semilogy(snr_vec, pbitDFE_same_channel)
semilogy(snr_vec_viterbi, pbit_viterbi_same_ch)
semilogy(snr_vec_fba_same_channel, pbit_fba_same_channel)
semilogy(snr_vec, BER_awgn(snr_vec))
semilogy(snr_vec, pbit_AWGN_sim(1:9))
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'Viterbi', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-7, 10^-1]), xlim([5, 15]), grid on
title('BER for various methods, same channel noise and estimate');


%% Statistics for estimated channel, average over independent experiments
snr_vec = 5:15;
load('pbit_viterbi')

figure
semilogy(snr_vec, mean(pbitLE, 2)), hold on,
semilogy(snr_vec, median(pbitDFE, 2))
semilogy(snr_vec_viterbi, mean(pbit_viterbi, 2))
semilogy(snr_vec, pbit_fba)
semilogy(snr_vec, BER_awgn(snr_vec))
semilogy(snr_vec, pbit_AWGN_sim)
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'Viterbi', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-5.5, 10^-1]), grid on
title('Average BER for independent experiments, estimated channel');
