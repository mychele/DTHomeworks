% Script to test the channel's polyphase implementation compared with the
% classic implementation

clear, clc, close all

Tc = 1;
T = 4*Tc;
snr = 1E30;

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
x = bitmap(trainingseq);    % Training symbols

Q0 = T/Tc; % interpolation factor
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];

% Output of the polyphase channel and of the classic implementation
out = channel_output(x, T, Tc, snr); % very high SNR!
out_control = conv(interp_nofilter(x, Q0), q);

% Make them the same length
difference = length(out) - length(out_control);
if difference > 0
    out_control = [out_control; zeros(difference, 1)];
elseif difference < 0
    out = [out; zeros(-difference, 1)];
end

% Plot comparison of the outputs
figure, stem(abs(out))
hold on, stem(abs(out_control), 'x')
title('Outputs of the channel compared')
legend('Polyphase output', 'Classic output')

% Compute and display error
err = sum(abs(out-out_control).^2) / sum(abs(out_control)).^2;
fprintf('The sum of the error squared, normalized wrt the energy of the\n')
fprintf('clean output, is %.2g = %.2f dB, with SNR=%gdB.\n', err, 10*log10(err), 10*log10(snr))