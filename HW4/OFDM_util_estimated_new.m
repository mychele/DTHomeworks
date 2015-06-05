%% ODFM script, it calls the OFDM function
clear
close all
clc
rng default

parpool(8);

%% Data
M = 512;
Npx = 7;
t0 = 5;
N2 = 4;
desired_bits = 2^22;
estMethod = 2;

%% BER with coding, estimated channel
repeat = 10;
isKnown = false;
coding = true;
snr_vec_coding_estimated_2 = [0, 1, 2:0.05:2.6, 3];
BER_coding_estimated_2 = zeros(length(snr_vec_coding_estimated_2), repeat);
parfor snr_i = 1:length(snr_vec_coding_estimated_2)
    snr_c = snr_vec_coding_estimated_2(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    for i = 1:repeat
        [BER_coding_estimated_2(snr_i, i), ~] = OFDM_BER(M, Npx, N2, t0, desired_bits, snr_c, coding, isKnown, estMethod);
    end
end

save('OFDM_coded_estimated_2', 'BER_coding_estimated_2', 'snr_vec_coding_estimated_2', 'desired_bits');
%% BER without coding, estimated channel
coding = false;
snr_vec_nocoding_estimated_2 = 0:15;
BER_nocoding_estimated_2 = zeros(length(snr_vec_nocoding_estimated_2), repeat);
parfor snr_i = 1:length(snr_vec_nocoding_estimated_2)
    snr_nc = snr_vec_nocoding_estimated_2(snr_i);
    fprintf('Uncoded, snr = %.2f\n', snr_nc);
    for i = 1:repeat 
        [BER_nocoding_estimated_2(snr_i, i), ~] = OFDM_BER(M, Npx, N2, t0, desired_bits, snr_nc, coding, isKnown, estMethod);
    end
end

save('OFDM_uncoded_estimated_2', 'BER_nocoding_estimated_2', 'snr_vec_nocoding_estimated_2', 'desired_bits');

delete(gcp);