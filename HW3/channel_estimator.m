% This script solves the first problem.
clear all
close all
clc

Tc = 1;
T = 4*Tc;

%% Training sequence generation
L = 15;
N = 10;
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
trainingseq = [mlseqdouble; mlseqdouble(1:2*N)];
trainingsymbols = bitmap(trainingseq);

%% Generate the channel output

snr = 20; %dB
snr_lin = 10^(snr/10);
r = channel_output(trainingsymbols, T, Tc, snr_lin);

%% Estimate qhat

