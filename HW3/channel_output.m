function [ r, sigma_w, q ] = channel_output( x, T, Tc, snr)
% CHANNEL_OUTPUT Generates channel output (that is the desired signal) via a
% polyphase implementation, with an hard coded non varying channel.
% Returns the channel output r given the input parameters
% x is a column vector for consistency
% snr must be linear

Q0 = T/Tc; % interpolation factor
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];
E_q = sum(abs(q).^2);
sigma_a = 2; % for a QPSK
sigma_w = sigma_a*E_q/(4*snr);

% create the impulse response for the Q0 = 4 branches
q_mat = [q(1:Q0:end); q(2:Q0:end); q(3:Q0:end); q(4:Q0:end)];

r = zeros(Q0 * (length(x)+length(q)/Q0-1), 1);

results = zeros(Q0, length(x)+length(q)/Q0-1);
for k = 1:Q0  % Iterate over the four branches
    w = wgn(1, length(x)+length(q)/Q0-1, 10*log10(sigma_w), 'complex');
    
    results(k,:) = conv(q_mat(k,:), x) + w.';
    r(k:Q0:end) = results(k,:);
end

end