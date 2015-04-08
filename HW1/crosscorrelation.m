function [ autoc ] = crosscorrelation( z1, z2, N_corr )
% CROSSCORRELATION - z1 and z2 should be of the same length
% This function is used only to find the crosscorrelation between complex
% sinuosoids in the frequency, amplitude and phase detector algorithm. Its
% definition recalls the definition of the xcorr function of MATLAB,
% although this is normalized.
K = length(z1);
autoc = zeros(N_corr + 1, 1);
for n = 1:(N_corr + 1)
    d = z1(n:K);
    b = conj(z2(1:(K - n + 1)));
    c = K - n + 1; 
    autoc(n) = d.' * b / c;
end

end