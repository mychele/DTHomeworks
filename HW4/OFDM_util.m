% ODFM transmitter

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

% Data
M = 512;
Npx = 7;
OFDM = true;
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

a = bitmap(int_enc_bits.');

%% Create data blocks

% perform a zero padding of the last symbols in order to have blocks of 512
% symbols
% TODO consider how can this affect BER
a_pad = [a; zeros(M - mod(length(a), M), 1)];
a_matrix = reshape(a_pad, M, []); % it should mantain columnwise order

% compute the ifft of blocks of 512 symbols
% http://it.mathworks.com/matlabcentral/newsreader/view_thread/17104
% it should be a columnwise operation!
A_matrix = ifft(a_matrix); % TODO check if matlab downscales
clear a_matrix % we care about memory consumption, don't we? 

% add the preamble to each column
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix]; % very powerful operation

% serialize in order to call channel output
s = reshape(A_matrix, [], 1);
% apply the correct gain too
s = M*s;

%% Send over the noisy channel
% because of 9.35, 9.38 and 9.39
snr = 14; % dB
snr_lin = 10^(snr/10);
[r, sigma_w, g] = channel_output(s, snr, OFDM);
G = fft(g, 512);
G = G(:);

%% Process at the receiver
r = r(t0:end); % consider the effect of the convolution at the end should be easy,
% since the resulting data should have a length which is a multiple of M +
% Npx
r = r(1:end - mod(length(r), M+Npx)); 

% perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);

G_i = G.^(-1);
x_matrix = fft(r_matrix);

y_matrix = bsxfun(@times, x_matrix, G_i);

y = reshape(y_matrix, [], 1);

% drop the zero padding
y = y(1:length(a));

%% Compute Log Likelihood Ratio
llr = zeros(2*length(y),1);
llr(1:2:end) = -2*real(y)/(sigma_w/2);
llr(2:2:end) = -2*imag(y)/(sigma_w/2);

%% Decode the bits

llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first

dec_bits = decodeBits(llr).';

Pbit = sum(xor(dec_bits, bits))/length(bits);








