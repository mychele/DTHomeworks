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
snr_vec_coding_known = 0.5:0.01:2.5;
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
%% BER with coding, estimated channel
repeat = 10;
isKnown = false;
coding = true;
snr_vec_coding_estimated = 2:0.1:3;
BER_coding_estimated = zeros(length(snr_vec_coding_estimated), repeat);
for snr_i = 1:length(snr_vec_coding_estimated)
    snr_c = snr_vec_coding_estimated(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    for i = 1:repeat
        [BER_coding_estimated(snr_i, i), ~] = OFDM_BER(M, Npx, desired_bits, snr_c, coding, isKnown);
    end
end

save('OFDM_coded_estimated', 'BER_coding_estimated', 'snr_vec_coding_estimated', 'desired_bits');
%% BER without coding, estimated channel
coding = false;
snr_vec_nocoding_estimated = 0:15;
BER_nocoding_estimated = zeros(length(snr_vec_nocoding_estimated), repeat);
for snr_i = 1:length(snr_vec_nocoding_estimated)
    snr_nc = snr_vec_nocoding_estimated(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    for i = 1:repeat 
        [BER_nocoding_estimated(snr_i, i), ~] = OFDM_BER(M, Npx, desired_bits, snr_nc, coding, isKnown);
    end
end

save('OFDM_uncoded_estimated', 'BER_nocoding_estimated', 'snr_vec_nocoding_estimated', 'desired_bits');
%% Plot
% figure,
% semilogy(snr_vec_coding_known, BER_coding),
% hold on,
% semilogy(snr_vec_nocoding_known, BER_nocoding)
