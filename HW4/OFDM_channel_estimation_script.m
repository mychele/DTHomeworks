%% Channel ESTIMATION for OFDM
% Send one block of data with symbols spaced of 16 channels

clear, close all
OFDM = true;
M = 512;
allowed_symb = 32;
spacing = M/allowed_symb;
Npx = 7;
t0 = 5;

block = ones(M, 1)*(-1-1i);

% Remember: for the symbols on which the estimation is performed, for a
% given snr (computed with the usual sigma_a^2 = 2), the power of the
% "estimation symbols" is doubled (-> better snr)
% Note that the variance of the noise at the receiver, after the DFT, is
% multplied by M, therefore it could be high
ts = ts_generation(allowed_symb-1, 1);
init_step = 1; % < 16
indices = init_step : spacing : init_step + spacing*(allowed_symb-1);
% Scale in order to double the power of tx symbols
% TODO What should we set the other symbols to?
block(indices) = ts * sqrt(2);

% Compute IDFT, add prefix, P/S
A = ifft(block);
A_pref = [A(end-Npx + 1:end); A];
s = reshape(A_pref, [], 1);

%% CHANNELIZATION

snr = 6; %dB
snr_lin = 10^(snr/10);
%fprintf('Symbols are pushed into the channel...\n');
% Send over the noisy channel
[r, sigma_w, g] = channel_output(s, snr_lin, OFDM);
g = g(1+t0 : end);   % Take t0 into account (just to plot stuff)
G = fft(g, 512);
G = G(:);


%% Process at the receiver

%fprintf('Symbols received, processing begins...\n');
r = r(1+t0 : end - mod(length(r), M+Npx) + t0);

% perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);
x_matrix = fft(r_matrix);

% Select useful samples
x_rcv = x_matrix(indices, 1);
X_known = diag(ts)*sqrt(2);

% LS estimaTION
phi = X_known'*X_known;
theta = X_known'*x_rcv;
G_est = phi \ theta;

% LS
F = dftmtx(M);
F = F(indices, 1:Npx+1);
% Solve LS for F*g=G_est where g is an 8x1 vector
g_hat = (F' * F) \ (F' * G_est);

g_est = ifft(G_est);
G_hat = fft(g_hat, M);


% Noise estimation
xhat = X_known * G_hat(indices);
E = sum(abs(xhat - x_rcv).^2)/length(xhat);
est_sigma_w = E/M;
fprintf('Est sigma_w^2 = %d\n', est_sigma_w);
fprintf('Real sigma_w^2 = %d\n', sigma_w);



% Plots

figure, hold on
stem(0:Npx, abs(g))
stem(0:Npx, abs(g_hat), 'x')
stem(0:15, abs(g_est(1:16)), '^')
legend('Actual g', 'g_hat', 'IDFT of G_est')

figure, 
subplot 211
plot(real(G)), hold on
plot(real(G_hat))
title(strcat('Comparison between estimated - LS+interpol - and real at ', num2str(snr), ' dB'))
legend('real(G)', 'real(G_hat)'), xlabel('i - subchannels'), ylabel('Real(G)'),
grid on, xlim([1, M])

subplot 212
plot(imag(G)), hold on
plot(imag(G_hat))
title(strcat('Comparison between estimated - LS+interpol - and real at ', num2str(snr), ' dB'))
legend('imag(G)', 'imag(G_hat)'), xlabel('i - subchannels'), ylabel('imag(G)'),
grid on, xlim([1, M])