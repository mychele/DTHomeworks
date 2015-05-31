function [packet, r, sigma_w] = txrc(L_data, snr)

% This script produces an output given by the 50 bit of ML training sequence and 
% L_data bit computed with an MLsequence (therefore L_data has to be 2^smth
% - 1)

%% Packet generation w/ ts, data
L = 31;
Nseq = 7;

trainingsymbols = ts_generation(L, Nseq);

MAX_ML = 2^20 - 1;
if (L_data > MAX_ML)
    dataseq = MLsequence(MAX_ML);
    dataseq = repmat(dataseq, ceil(L_data / MAX_ML), 1);
    dataseq = dataseq(1:L_data);
else
    dataseq = MLsequence(L_data);
end
datasymbols = bitmap(dataseq(1 : end - mod(L_data, 2)));
% QPSK requires an even number of bits

packet = [trainingsymbols; datasymbols];
%% Generate the channel output
snr_lin = 10^(snr/10);
[r, sigma_w, ~] = channel_output(packet, snr_lin, false);

end