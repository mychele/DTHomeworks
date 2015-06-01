%% Channel ESTIMATION for OFDM
% Very rough script
% Send one block of data with symbols spaced of 16 channels

OFDM = true;
M = 512;
allowed_symb = 32;
spacing = M/allowed_symb;
Npx = 7;

symbol = 1+1i;
%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%
% NOT SURE IF IT'S OK TO SET TO ZERO THE OTHER SYMBOLS!!
block = zeros(M, 1);
block(1:spacing:end) = symbol;

A = ifft(block);
A_pref = [A(end-Npx + 1:end); A];

%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%
% NOT SURE IF IT'S OK TO SCALE IN THIS WAY
% Remember: for the symbols on which the estimation is performed, for a
% given snr (computed with the usual sigma_a^2 = 2), the power of the
% "estimation symbols" is doubled (-> better snr)
% Note that the variance of the noise at the receiver, after the DFT, is
% multplied by M, therefore it could be high
% Scale in order to double the power of tx symbols
A_pref = A_pref*sqrt(2);

s = reshape(A, [], 1);

%% CHANNELIZATION
snr = 100; %dB
snr_lin = 10^(snr/10);
fprintf('Symbols are pushed into the channel...\n');
% Send over the noisy channel
[r, sigma_w, g] = channel_output(s, snr_lin, OFDM);
G = fft(g, 512);
G = G(:);

%% Process at the receiver
fprintf('Symbols received, processing begins...\n');
r = r(1:end - mod(length(r), M+Npx));

% perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);
x_matrix = fft(r_matrix);

% LS estimaTION
% Select useful samples
x_rcv = x_matrix(1:spacing:end, 1);
X_known = diag(symbol*ones(length(x_rcv), 1));

phi = X_known'*X_known;
theta = X_known'*x_rcv;

G_est = phi \ theta;

% InterpolaTION
f = 1:spacing:M;
f_fine = 1:1:M;
G_est_complete = interp1(f, G_est, f_fine, 'pchip');

figure, hold on
plot(20*log10(abs(G_est_complete))), 
plot(1:spacing:M, 20*log10(abs(G_est)), '-h'),
plot(20*log10(abs(G))),
title(strcat('Comparison between estimated - LS+interpol - and real at ', num2str(snr), ' dB'))
legend('|G_{est}| iterpolated', '|G_{est}| not interpolated', '|G|'), xlabel('i - subchannels'), ylabel('|G| - dB'),
grid on



