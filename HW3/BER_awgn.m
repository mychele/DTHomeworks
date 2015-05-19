function [BER] = BER_awgn(snr)


snr_lin = 10.^(snr./10);

%% BER for ideal AWGN channel

q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];

h = q(3:4:end);

E_h = sum(abs(h).^2);
dm = 2; % for a QPSK, it is also the distance between 2 symbols
sigma_w_2 = 2*E_h./(snr_lin);
sigma_i = sqrt(sigma_w_2/2);

gamma = (dm ./ (2*sigma_i)).^2;

BER = qfunc(sqrt(gamma));

end