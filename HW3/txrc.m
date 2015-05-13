function [r, sigma_w] = txrc(L_data, snr, T, Tc)

% This script produces an output given by the 50 bit of ML training sequence and 
% L_data bit computed with an MLsequence (therefore L_data has to be 2^smth
% - 1)

%% Training sequence generation
L = 15;
Nseq = 10;
mlseq = MLsequence(L);

% Replace every 0 with two 0s and every 1 with two 1s to put it into the
% bitmap
mlseqdouble = zeros(2*L,1);
for i = 1:L
    switch mlseq(i)
        case 0
            mlseqdouble(2*i-1) = 0;
            mlseqdouble(2*i) = 0;
        case 1
            mlseqdouble(2*i-1) = 1;
            mlseqdouble(2*i) = 1;
    end
end

% Repeat the sequence and bitmap it to get the symbols
trainingseq = [mlseqdouble; mlseqdouble(1:2*Nseq)];
trainingsymbols = bitmap(trainingseq);

dataseq = MLsequence(L_data);
datasymbols = bitmap(dataseq(1:end-1)); % the ML sequence has an odd length,
% qpsk requires even symbols

packet = [trainingsymbols; datasymbols];

%% Generate the channel output
snr_lin = 10^(snr/10);
[r, sigma_w, ~] = channel_output(packet, T, Tc, snr_lin);

end