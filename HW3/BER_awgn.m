function [BER] = BER_awgn(snr)

snr_lin = 10.^(snr./10);
BER = qfunc(sqrt(snr_lin));

end