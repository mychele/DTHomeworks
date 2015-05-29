% This script converts bits into symbols, using the encoder + interleaver +
% bitmap setup.
clear all
close all
clc

rng default

%% Get optimal number of bits
desired_bits = 2^18;
% Compute the closest number of bits that both interleaver and encoder will
% like
found = false;
bit_number = 0;
while(~found)
    search_step = 32400;
    bit_number = bit_number + search_step;
    if (bit_number > desired_bits)
        found = true;
    end
end

%% Generate and encode bits
bits = randi([0 1], 1, bit_number);

enc_bits = encodeBits(bits);

int_enc_bits = interleaver(enc_bits);  % Interleave the encoded bits

symbols = bitmap(int_enc_bits.'); 

%% Send stuff through
snrdb = 2.3;
snrlin = 10^(snrdb/10);
[rcv_bits, sigma_w, h] = channel_output(symbols, snrlin);

rcv_bits = rcv_bits(6:end-7);

% Eh = 1;

% sigma_w = 2*Eh/snrlin;
% w = wgn(length(symbols), 1, 10*log10(sigma_w), 'complex');
% rcv_bits = symbols + w;

%% Receiver: filter with DFE
[Jmin, rcv_bits] = DFE_filter(rcv_bits, h.', 5, 7, sigma_w, 15, 20, 7+20-1-15, false);

%% Compute Log Likelihood Ratio
llr = zeros(2*length(symbols),1);
llr(1:2:end) = -2*real(rcv_bits)/(sigma_w/2);
llr(2:2:end) = -2*imag(rcv_bits)/(sigma_w/2);

%% Decode the bits

llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first

dec_bits = decodeBits(llr).';

Pbit = sum(xor(dec_bits, bits))/length(bits);