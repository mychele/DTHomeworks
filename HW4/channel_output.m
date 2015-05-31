function [ r, sigma_w, h ] = channel_output( x, snr, OFDM)
% CHANNEL_OUTPUT Generates channel output (that is the desired signal)
% with an hard coded non varying channel.
% Returns the channel output r given the input parameters
% x is a column vector for consistency
% snr must be linear

h = [0,0,0,0,0,0.7*exp(-1i*2.57), 0.24*exp(1i*1.34), 0.15*exp(-1i*2.66), ...
    0.58*exp(-1i*1.51), 0.4*exp(-1i*1.63), 0, 0, 0];
M = 512;

E_h = sum(abs(h).^2);
if (OFDM == true) % different definition of SNR Msigma_a E_h / sigma_w
    sigma_a = 2/M;
else
    sigma_a = 2; % for a QPSK
end
sigma_w = sigma_a*E_h/snr;

w = wgn(1, length(x) + length(h) - 1, 10*log10(sigma_w), 'complex');

r = conv(h, x) + w.';
end