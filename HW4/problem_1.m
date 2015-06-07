% This script solves the first problem.
clear
close all
clc

rng default

%parpool(15);

% data
sigma_a = 2;

%% Get optimal number of bits
desired_bits = 2^24;
% Compute the closest number of bits that both interleaver and encoder will like
search_step = 32400;
bit_number = ceil(desired_bits / search_step) * search_step;

%% Estimate Pbit for ideal channel without encoding
snr_vec = 0:0.5:14;
Pbit_noenc = zeros(1,length(snr_vec));

parfor curr_snr = 1:length(snr_vec)
   disp(curr_snr);
   snrdb = snr_vec(curr_snr);
   % Simulate Pbit for the ideal channel, no encoding
   bits = randi([0 1], 1, bit_number);
   symbols = bitmap(bits.');
   
   % Send data through the ideal channel
   snrlin = 10^(snrdb/10);
   Eh = 1;  % Energy of the ideal channel ir
   sigma_w = sigma_a*Eh/snrlin;
   w = wgn(length(symbols), 1, 10*log10(sigma_w), 'complex');
   rcv_bits = symbols + w;
   
   % Threshold the bits
   decided_symbols = zeros(1, length(rcv_bits));
   for idx = 1:length(rcv_bits)
      decided_symbols(idx) = qpsk_td(rcv_bits(idx));
   end
   
   decided_bits = ibmap(decided_symbols);
   
   % Estimate pbit
   Pbit_noenc(curr_snr) = sum(xor(decided_bits.', bits))/length(bits);
   
end
   
%% Estimate Pbit for ideal channel with encoding
snr_vec_enc = 0:0.02:1;
Pbit_enc = zeros(1,length(snr_vec_enc));

parfor curr_snr = 1:length(snr_vec_enc)
   disp(snr_vec_enc(curr_snr));
   snrdb = snr_vec_enc(curr_snr);
   bits = randi([0 1], 1, bit_number);
   
   enc_bits = encodeBits(bits);
   int_enc_bits = interleaver(enc_bits);  % Interleave the encoded bits
   
   symbols = bitmap(int_enc_bits.');
   
   % Send the data through the ideal channel
   snrlin = 10^(snrdb/10);
   Eh = 1;
   sigma_w = sigma_a*Eh/snrlin;
   w = wgn(length(symbols), 1, 10*log10(sigma_w), 'complex');
   rcv_bits = symbols + w;
   
   % Compute Log Likelihood Ratio
   llr = zeros(2*length(symbols),1);
   llr(1:2:end) = -2*real(rcv_bits)/(sigma_w/2);
   llr(2:2:end) = -2*imag(rcv_bits)/(sigma_w/2);
   
   % Decode the bits
   llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first
   dec_bits = decodeBits(llr).';
   
   Pbit_enc(curr_snr) = sum(xor(dec_bits, bits))/length(bits);
   
end

%% Save results
save('Problem1', 'snr_vec', 'snr_vec_enc', 'bit_number', 'Pbit_noenc', 'Pbit_enc');

%% Plot results
load ('Problem1');
semilogy(snr_vec, Pbit_noenc), hold on,
semilogy(snr_vec_enc, Pbit_enc)
ylim([10^-5, 10^-1]), xlim([0, 14]), grid on
xlabel('\Gamma (dB)')
ylabel('Pbit')
legend('Uncoded QPSK', 'Coded QPSK')
title('Bit Error Rate for uncoded and coded QPSK');