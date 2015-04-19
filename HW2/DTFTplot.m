function [ ] = DTFTplot(b, a, res , fs)
%DTFTPLOTLOG Logarithmic plot of the DTFT module of the signal or filter
%with numerator b and denominator a, res is the sampling freq

[Hf, f] = freqz(b, a, res, 'whole', fs);

figure
plot(f, 20*log10(abs(Hf)))

end