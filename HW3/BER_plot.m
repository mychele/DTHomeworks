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
semilogy(snr_vec, pbitDFE(:, 2)), hold on, 
semilogy(snr_vec_viterbi, pbit_viterbi(:, 3)), hold on,
semilogy(snr_vec, pbit_fba), hold on,
semilogy(snr_vec, BER_awgn(snr_vec)), hold on,
semilogy(snr_vec, pbit_AWGN_sim)
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'Viterbi', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-5.5, 10^-1]), grid on


%% Statistics for real channel

load('pbit_LE_real_channel.mat');
load('pbit_DFE_real_channel.mat');
load('pbit_fba_real_channel.mat');
load('pbit_viterbi_realch.mat')


figure
semilogy(snr_vec, pbitLE_real_channel), hold on,
semilogy(snr_vec, pbitDFE_real_channel), hold on, 
semilogy(snr_vec_viterbi + 0.3, pbit_viterbi), hold on,
semilogy(snr_vec_fba_real_channel + 0.3, pbit_fba_real_channel), hold on,
semilogy(snr_vec, BER_awgn(snr_vec)), hold on,
semilogy(snr_vec, pbit_AWGN_sim)
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-7, 10^-1]), xlim([5, 15]), grid on

%% Statistics for same estimate

snr_vec = 5:13; 

load('pbit_LE_same_channel.mat');
load('pbit_DFE_same_channel.mat');
load('pbit_viterbi_same_ch513.mat');

figure
semilogy(snr_vec, pbitLE_same_channel), hold on,
semilogy(snr_vec, pbitDFE_same_channel), hold on, 
semilogy(snr_vec_viterbi, pbit_viterbi_same_ch), hold on,
%semilogy(snr_vec_fba_real_channel + 0.3, pbit_fba_real_channel), hold on,
semilogy(snr_vec, BER_awgn(snr_vec)), hold on,
semilogy(snr_vec, pbit_AWGN_sim(1:9))
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 15', 'DFE, M1 = 25, D = 24, M2 = 4', 'FBA', 'AWGN', 'AWGN with simulation')
ylim([10^-7, 10^-1]), xlim([5, 15]), grid on