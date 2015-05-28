% This script converts bits into symbols, using the encoder + interleaver +
% bitmap setup.
clear all
close all
clc

rng default

%% Initialization
% Encoder
enc = fec.ldpcenc;

% Decoder
dec = fec.ldpcdec;
dec.DecisionType = 'Hard Decision';
dec.OutputFormat = 'Information Part';
dec.NumIterations = 50;
dec.DoParityChecks = 'Yes';

%% Generate bits using an ML sequence
% L = 31;
% bits = MLsequence(L).';

%% Encode the bits
bits = randi([0 1], 1, enc.NumInfoBits); % The encoder wants NumInfoBits bits

enc_bits = encode(enc, bits);   % length(enc_bits) = 2*length(bits)

int_enc_bits = interleaver(enc_bits);  % Interleave the encoded bits

symbols = bitmap(int_enc_bits.'); 

%% Send stuff through
Eh = 1;
sigma_w = 0.0001;
w = wgn(length(symbols), 1, 10*log10(sigma_w), 'complex');
rcv_bits = symbols + w;

%% Compute Log Likelihood Ratio
llr = zeros(64800,1);
llr(1:2:end) = -2*real(rcv_bits)/(sigma_w/2);
llr(2:2:end) = -2*imag(rcv_bits)/(sigma_w/2);

%% Decode the bits

llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first

dec_bits = decode(dec, llr);

Pbit = sum(xor(dec_bits, bits))/length(bits);