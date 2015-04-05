function [ ] = DTFTplot( x, res )
%DTFTPLOTLOG Summary of this function goes here
%   Detailed explanation goes here

b = x; a = 1;  % Define b and a in z^-1
[Hf, f] = freqz(b, a, res, 1, 'whole');

figure
%subplot(2, 1, 1)
plot(f, 10*log10(abs(Hf)))
%subplot(2, 1, 2), plot(f, angle(Hf))

end