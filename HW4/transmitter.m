% This script converts bits into symbols, using the encoder + interleaver +
% bitmap setup.
clear
close all
clc

rng default
% The design of the DFE equalizer has to be carried out assuming the
% channel is known
% From the assigned impulse response
t0 = 6; 
N1 = 0;
N2 = 4;
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
coding = true;

symbols = bitmap(int_enc_bits.'); 

%% Send stuff through
snrdb = 2.15;
snrlin = 10^(snrdb/10);
[rcv_symb, sigma_w, h] = channel_output(symbols, snrlin);

% Normalization!
rcv_symb = rcv_symb(t0:end-7)/h(t0);
hi = h(t0-N1:t0+N2)/h(t0);

%% Receiver: filter with DFE
M1_dfe = 15;
D_dfe = M1_dfe - 1;
M2_dfe = N2 + M1_dfe - 1 - D_dfe;
[~, rcv_bits] = DFE_filter(rcv_symb, hi.', N1, N2, sigma_w, D_dfe, M1_dfe, M2_dfe, coding, false);

%% Compute Log Likelihood Ratio
llr = zeros(2*length(symbols),1);
llr(1:2:end) = -2*real(rcv_bits)/(sigma_w/2);
llr(2:2:end) = -2*imag(rcv_bits)/(sigma_w/2);

%% Decode the bits

llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first

dec_bits = decodeBits(llr).';

Pbit = sum(xor(dec_bits, bits))/length(bits);