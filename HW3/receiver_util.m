% This script produces data that can be used by the receiver to perform
% detection. Use this to simulate tx with different SNRs.
clear
close all
clc
% rng default

T = 1;
snr = 14; % 6, 8, 10, 12, 14 % dB
L_data = 2^20-1; %2^18-1;

% From exercise 1
assumed_dly = 2;
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);  % offset in T/4
t0 = assumed_dly;                    % t0 is @ T; TODO this is N1, we should refactor

%% Create, send and receive data, estimate channel and prepare for detection

% Create, send and receive data with the given channel
fprintf('Generating input symbols and channel output... ')
[packet, r, sigma_w] = txrc(L_data, snr, assumed_m_opt);
fprintf('done!\n')

% Estimate the channel using the first 100 samples (4*length(ts))
fprintf('Estimating timing phase and IR... ')
[ h, est_sigmaw, N1, N2 ] = get_channel_info(r(init_offs+1:25+init_offs), 0, 4, assumed_m_opt);
fprintf('done!\n')

% Sample to get r @ T
x = r / h(N1+1).';         % data normalized by h0
hi = h / h(N1+1).';         % impulse response normalized by h0

%% Detection begins

N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1 = 20;   % FF filter: equal to the span of h
D = 14; %M1-1;   % D is chosen large first and then decreased % (N-1)/2 + 2 or 3 for LE
M2 = 0; %N2 + M1 - 1 - D;      % FB filter: one less than the FF filter
verb = 1;
if (L_data > 128) % to avoid very looooooong time to generate useless plots
    verb=0;
end
[decisions, pbit, num_bit_error, Jmin] = DFE_filter(packet, x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2, est_sigmaw, t0, D, M1, M2, verb);
%[decisions, pbit, num_bit_error] = viterbi(packet(1:end-assumed_dly), x(1+assumed_dly:end), hi, N1, N2, 0, N2);


