function [ BER, G ] = OFDM_BER( M, Npx, t0, desired_bits, snr, coding, chIsKnown, varargin )
%This function performs the transmission and reception of bits with ODFM,
%with or without encoding
%   It needs
%   - M block length
%   - Npx prefix length
%   - the number of desired bits
%   - the snr
%   - the option coding to be set either true or false
%   - the estimation method (1 standard, 2 8 points+noise) - optional, if
%   not indicated the classic is used
%   It returns
%   - the BER
%   - the channel frequency response
%   This implementation doesn't rescale the power of sent symbols after the
%   IFFT operation, therefore sigma_s^2 = sigma_a^/M

if (length(varargin) == 1)
    if (varargin{1} == 1)
        fprintf('Use classic estimation method \n')
        estMethod = 1;
    elseif (varargin{1} == 2)
        fprintf('Use new estimation method \n')
        estMethod = 2;
    end
else
    fprintf('Use classic estimation method \n')
    estMethod = 1;
end

warning('off', 'all');
OFDM = true;

% Compute the optimal number of bits
fprintf('Start transmission...\n');
if (coding == true)
    % Compute the closest number of bits that both interleaver and encoder will like
    search_step = 32400;
    bit_number = ceil(desired_bits / search_step) * search_step;
else
    bit_number = desired_bits;
end

% Generate and encode bits
bits = randi([0 1], 1, bit_number).';

if (coding == true)
    enc_bits = encodeBits(bits.');
    int_enc_bits = interleaver(enc_bits);  % Interleave the encoded bits
    a = bitmap(int_enc_bits.');
else
    a = bitmap(bits);
end

% Create data blocks

% perform a padding of the last symbols in order to have blocks of 512 symbols: we use
% -1-j to perform the padding
a_pad = [a; ones(M - mod(length(a), M), 1) * (-1-1i)];
a_matrix = reshape(a_pad, M, []); % it should mantain columnwise order

% compute the ifft of blocks of 512 symbols
% http://it.mathworks.com/matlabcentral/newsreader/view_thread/17104
% it should be a columnwise operation!
A_matrix = ifft(a_matrix);

% add the preamble to each column
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix]; % very powerful operation

% serialize in order to call channel output
s = reshape(A_matrix, [], 1);

fprintf('Symbols are pushed into the channel...\n');
% Send over the noisy channel
snr_lin = 10^(snr/10);
[r, sigma_w, g] = channel_output(s, snr_lin, OFDM);
if (chIsKnown)
    g = g(1+t0 : end);
    G = fft(g, 512);
    G = G(:);
else
    if(estMethod == 1)
        [G, sigma_w] = OFDM_channel_estimation(snr, Npx, t0);
    elseif(estMethod == 2)
        [G, sigma_w] = OFDM_channel_estimation_2(snr, Npx, t0);   
    end 
end

% Process at the receiver
% consider the effect of the convolution at the end should be easy,
% since the resulting data should have a length which is a multiple of M +
% Npx
fprintf('Symbols received, processing begins...\n');
r = r(1+t0 : end - mod(length(r), M+Npx) + t0);

% perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);

G_inv = G.^(-1);
x_matrix = fft(r_matrix);

y_matrix = bsxfun(@times, x_matrix, G_inv);

% Detect and compute BER
if (coding == true)
    % sigma_i after the DFT and the scaling by G_i of each branch
    sigma_i = 0.5*sigma_w*M*abs(G_inv).^2;
    % Compute Log Likelihood Ratio
    % It is different for each branch
    llr_real = -2*bsxfun(@times, real(y_matrix), sigma_i.^(-1));
    llr_imag = -2*bsxfun(@times, imag(y_matrix), sigma_i.^(-1));
    llr_real_ar = reshape(llr_real, [], 1);
    llr_imag_ar = reshape(llr_imag, [], 1);
    llr = zeros(numel(llr_real) + numel(llr_imag), 1);
    llr(1:2:end) = llr_real_ar;
    llr(2:2:end) = llr_imag_ar;
    % Drop the zero padding
    llr = llr(1:length(enc_bits));
    % Decode the bits
    llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first
    dec_bits = decodeBits(llr).';
else
    y = reshape(y_matrix, [], 1);
    % drop the zero padding
    y = y(1:length(a));
    decision = zeros(length(y), 1);
    for k = 1:length(y)
        decision(k) = qpsk_td(y(k));
    end
    dec_bits = ibmap(decision);
end

% make dec_bits a column even if it is already
dec_bits = dec_bits(:);

BER = sum(xor(dec_bits, bits))/length(bits);
fprintf('End, the BER is %d \n', BER);

% my_llr_linear = 10.^(llr(1:100)/10);
% my_dec_bits = dec_bits(1:100);
% my_p = 1 ./ (1 + my_llr_linear);
% figure, stem(my_p), hold on, stem(my_dec_bits)

end
