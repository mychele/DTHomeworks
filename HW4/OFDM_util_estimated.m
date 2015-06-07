%% ODFM script, it calls the OFDM function
clear
close all
clc
rng default

%parpool(8);

%% Data
M = 512;
Npx = 7;
t0 = 5;
N2 = 4;
desired_bits = 2^22;

%% BER with coding, estimated channel
repeat = 50;
isKnown = false;
coding = true;
snr_vec_coding_estimated_1 = 2:0.05:2.6;
BER_coding_estimated_1 = zeros(length(snr_vec_coding_estimated_1), repeat);
for snr_i = 1:length(snr_vec_coding_estimated_1)
    snr_c = snr_vec_coding_estimated_1(snr_i);
    fprintf('Coded, snr = %.2f\n', snr_c);
    for i = 1:repeat
        [BER_coding_estimated_1(snr_i, i), ~] = OFDM_BER(M, Npx, N2, t0, desired_bits, snr_c, coding, isKnown);
    end
end

save('OFDM_coded_estimated_long_1', 'BER_coding_estimated_1', 'snr_vec_coding_estimated_1', 'desired_bits');
%% BER without coding, estimated channel
% coding = false;
% snr_vec_nocoding_estimated = 0:15;
% BER_nocoding_estimated = zeros(length(snr_vec_nocoding_estimated), repeat);
% parfor snr_i = 1:length(snr_vec_nocoding_estimated)
%     snr_nc = snr_vec_nocoding_estimated(snr_i);
%     fprintf('Uncoded, snr = %.2f\n', snr_nc);
%     for i = 1:repeat 
%         [BER_nocoding_estimated(snr_i, i), ~] = OFDM_BER(M, Npx, N2, t0, desired_bits, snr_nc, coding, isKnown);
%     end
% end
% 
% save('OFDM_uncoded_estimated_long', 'BER_nocoding_estimated', 'snr_vec_nocoding_estimated', 'desired_bits');
% 
% delete(gcp);