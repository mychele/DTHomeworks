%% ODFM script, it calls the OFDM function and performs different BER estimates
clear
close all
clc
rng default

%% Data
M = 512;
Npx = 7;
desired_bits = 2^22;
estMethod = 2;
t0 = 5;

%% BER with coding, known channel

isKnown = true;
coding = true;
snr_vec_coding_known = 0:0.05:4; % no more is needed
BER_coding_known = zeros(length(snr_vec_coding_known), 1);
for snr_i = 1:length(snr_vec_coding_known)
    snr_c = snr_vec_coding_known(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    [BER_coding_known(snr_i), ~] = OFDM_BER(M, Npx, t0, desired_bits, snr_c, coding, isKnown);
end

save('OFDM_coded_known', 'BER_coding_known', 'snr_vec_coding_known', 'desired_bits');

%% BER without coding, known channel

coding = false;
snr_vec_nocoding_known = 0:15;
BER_nocoding_known = zeros(length(snr_vec_nocoding_known), 1);
for snr_i = 1:length(snr_vec_nocoding_known)
    snr_nc = snr_vec_nocoding_known(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    [BER_nocoding_known(snr_i), ~] = OFDM_BER(M, Npx, t0, desired_bits, snr_nc, coding, isKnown);
end

save('OFDM_uncoded_known', 'BER_nocoding_known', 'snr_vec_nocoding_known', 'desired_bits');

%% BER with coding, estimated channel

isKnown = false;
coding = true;
snr_vec_coding_estimated = [0, 1, 2:0.05:3];
BER_coding_estimated = zeros(length(snr_vec_coding_estimated), 1);
for snr_i = 1:length(snr_vec_coding_estimated)
    snr_c = snr_vec_coding_estimated(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    [BER_coding_estimated(snr_i), ~] = OFDM_BER(M, Npx, t0, desired_bits, snr_c, coding, isKnown, estMethod);
    
end

save('OFDM_coded_estimated', 'BER_coding_estimated', 'snr_vec_coding_estimated', 'desired_bits');

%% BER without coding, estimated channel

coding = false;
snr_vec_nocoding_estimated = 0:15;
BER_nocoding_estimated = zeros(length(snr_vec_nocoding_estimated), 1);
for snr_i = 1:length(snr_vec_nocoding_estimated)
    snr_nc = snr_vec_nocoding_estimated(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    [BER_nocoding_estimated(snr_i), ~] = OFDM_BER(M, Npx, t0, desired_bits, snr_nc, coding, isKnown, estMethod); 
end

save('OFDM_uncoded_estimated', 'BER_nocoding_estimated', 'snr_vec_nocoding_estimated', 'desired_bits');

%% BER with coding, estimated channel, first method

estMethod = 1;
isKnown = false;
coding = true;
snr_vec_coding_estimated = [0, 1, 2:0.05:3];
BER_coding_estimated = zeros(length(snr_vec_coding_estimated), 1);
for snr_i = 1:length(snr_vec_coding_estimated)
    snr_c = snr_vec_coding_estimated(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    [BER_coding_estimated(snr_i), ~] = OFDM_BER(M, Npx, t0, desired_bits, snr_c, coding, isKnown, estMethod);
    
end

save('OFDM_coded_estimated_1', 'BER_coding_estimated', 'snr_vec_coding_estimated', 'desired_bits');

%% BER without coding, estimated channel, first method

coding = false;
snr_vec_nocoding_estimated = 0:15;
BER_nocoding_estimated = zeros(length(snr_vec_nocoding_estimated), 1);
for snr_i = 1:length(snr_vec_nocoding_estimated)
    snr_nc = snr_vec_nocoding_estimated(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    [BER_nocoding_estimated(snr_i), ~] = OFDM_BER(M, Npx, t0, desired_bits, snr_nc, coding, isKnown, estMethod); 
end

save('OFDM_uncoded_estimated_1', 'BER_nocoding_estimated', 'snr_vec_nocoding_estimated', 'desired_bits');
