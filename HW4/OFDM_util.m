%% ODFM util script, it calls the OFDM function
clear
close all
clc
rng default

parpool(4); 

%% Data
M = 512;
Npx = 7;
desired_bits = 2^20;

%% BER with coding
coding = true;
snr_vec_coding = 1:0.05:2;
BER_coding = zeros(length(snr_vec_coding), 1);
parfor snr_i = 1:length(snr_vec_coding)
    snr_c = snr_vec_coding(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    [BER_coding(snr_i), ~] = OFDM_BER(M, Npx, desired_bits, snr_c, coding);
end

%% BER without coding
coding = false;
snr_vec_nocoding = 0:15;
BER_nocoding = zeros(length(snr_vec_nocoding), 1);
parfor snr_i = 1:length(snr_vec_nocoding)
    snr_nc = snr_vec_nocoding(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    [BER_nocoding(snr_i), ~] = OFDM_BER(M, Npx, desired_bits, snr_nc, coding);
end

delete(gcp);
%% Plot

figure, 
semilogy(snr_vec_coding, BER_coding),
hold on,
semilogy(snr_vec_nocoding, BER_nocoding)












