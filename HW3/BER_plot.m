%% Group all the BER curves

load('pbit_DFE.mat');
load('pbit_LE.mat');
load('pbit_viterbi.mat');

%% Statistics
snr_vec = 6:2:14;

BER_LE = mean(pbitLE, 2);
BER_DFE = mean(pbitDFE, 2);
BER_viterbi = mean(pbit_viterbi, 2);

BER_ideal = BER_awgn(snr_vec);

figure, hold on
semilogy(snr_vec, BER_LE)
semilogy(snr_vec, BER_DFE)
semilogy(snr_vec_viterbi, BER_viterbi)
semilogy(snr_vec, BER_ideal)
xlabel('snr [dB]'), ylabel('BER')
legend('LE, M1 = 20, D = 16', 'DFE, M1 = 25, D = 24', 'Viterbi', 'AWGN')
ylim([10^-6, 10^-1])

figure, semilogy(snr_vec, BER_ideal),
xlabel('snr [dB]'), ylabel('BER'), legend( 'AWGN')
ylim([10^-6, 10^-1])