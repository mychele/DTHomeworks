function [ r ] = channel_output( x, T, Tc, sigma_w)
% CHANNEL_OUTPUT Generates channel output (that is the desired signal) via a
% polyphase implementation, with an hard coded non varying channel.
% Returns the channel output r given the input parameters
% x is a row vector

Q0 = T/Tc; % interpolation factor
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];

% create the impulse response for the Q0 = 4 branches
q_mat = [q(1:Q0:end); q(2:Q0:end); q(3:Q0:end); q(4:Q0:end)];

r = zeros(T * length(x), 1);

for k = 0 : length(x)-1
    % Generate white noise for each branch
    w = wgn(Q0, 1, 10*log10(sigma_w), 'complex');
    if (k < length(q)/Q0)
        xconv = [fliplr(x(1:k+1)), zeros(1, length(q)/Q0 - k - 1)];
    else
        xconv = fliplr(x(k-length(q)/Q0 + 2:k+1));
    end
    for i = 0:Q0-1 % Branch index
        r(k*T + i*Tc + 1) = q_mat(i + 1, :) * xconv.' + w(i+1);
    end
end

end