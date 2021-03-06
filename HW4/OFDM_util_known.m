%% ODFM script, it calls the OFDM function
clear
close all
clc
rng default

%% Data
M = 512;
Npx = 7;
desired_bits = 2^22;

%% BER with coding, known channel
isKnown = true;
coding = true;
snr_vec_coding_known = 0.5:0.05:2;
BER_coding_known = zeros(length(snr_vec_coding_known), 1);
for snr_i = 1:length(snr_vec_coding_known)
    snr_c = snr_vec_coding_known(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    [BER_coding_known(snr_i), ~] = OFDM_BER(M, Npx, desired_bits, snr_c, coding, isKnown);
end

save('OFDM_coded_known', 'BER_coding_known', 'snr_vec_coding_known', 'desired_bits');
%% BER without coding, known channel
coding = false;
snr_vec_nocoding_known = 0:15;
BER_nocoding_known = zeros(length(snr_vec_nocoding_known), 1);
for snr_i = 1:length(snr_vec_nocoding_known)
    snr_nc = snr_vec_nocoding_known(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    [BER_nocoding_known(snr_i), ~] = OFDM_BER(M, Npx, desired_bits, snr_nc, coding, isKnown);
end

save('OFDM_uncoded_known', 'BER_nocoding_known', 'snr_vec_nocoding_known', 'desired_bits');

%% Plot
% figure,
% semilogy(snr_vec_coding_known, BER_coding),
% hold on,
% semilogy(snr_vec_nocoding_known, BER_nocoding)
