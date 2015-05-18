function [ r, sigma_w, h ] = channel_output_T( x, snr, m_opt)
% CHANNEL_OUTPUT Generates channel output (that is the desired signal)
% with an hard coded non varying channel.
% Returns the channel output r given the input parameters
% x is a column vector for consistency
% snr must be linear
% m_opt is the optimal delay in T/4

Q0 = 4; % interpolation factor
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];
init_offs = mod(m_opt, 4);
h = q(init_offs + 1:Q0:end);

E_h = sum(abs(h).^2);
sigma_a = 2; % for a QPSK
sigma_w = sigma_a*E_h/snr;

w = wgn(1, length(x) + length(h) - 1, 10*log10(sigma_w), 'complex');

r = conv(h, x) + w.';
end