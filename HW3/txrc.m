function [packet, r, sigma_w] = txrc(L_data, snr, m_opt)

% This script produces an output given by the 50 bit of ML training sequence and 
% L_data bit computed with an MLsequence (therefore L_data has to be 2^smth
% - 1)

%% Packet generation w/ ts, data
L = 15;
Nseq = 10;

trainingsymbols = ts_generation(L, Nseq);

dataseq = MLsequence(L_data);
datasymbols = bitmap(dataseq(1:end-1)); % the ML sequence has an odd length,
% qpsk requires even symbols

packet = [trainingsymbols; datasymbols];
%% Generate the channel output
snr_lin = 10^(snr/10);
[r, sigma_w, ~] = channel_output_T(packet, snr_lin, m_opt);

end