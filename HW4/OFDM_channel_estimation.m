function [G_hat, est_sigma_w] = OFDM_channel_estimation(snr, Npx, N2, t0)

%% Channel ESTIMATION for OFDM
% Send one block of data with symbols spaced of 16 channels

OFDM = true;
M = 512;
allowed_symb = 32;
spacing = M/allowed_symb;

block = ones(M, 1)*(-1-1i);

% Remember: for the symbols on which the estimation is performed, for a
% given snr (computed with the usual sigma_a^2 = 2), the power of the
% "estimation symbols" is doubled (-> better snr)
% Note that the variance of the noise at the receiver, after the DFT, is
% multplied by M, therefore it could be high
ts = ts_generation(allowed_symb-1, 1) * sqrt(2);
init_step = 1; % < 16
indices = init_step : spacing : init_step + spacing*(allowed_symb-1);
% Scale in order to double the power of tx symbols
% TODO What should we set the other symbols to?
block(indices) = ts;

% Compute IDFT, add prefix, P/S
A = ifft(block);
A_pref = [A(end-Npx + 1:end); A];
s = reshape(A_pref, [], 1);

%% CHANNELIZATION
% Do our own transmission and reception of the training sequence

snr_lin = 10^(snr/10);
[r, ~, ~] = channel_output(s, snr_lin, OFDM);


%% Process at the receiver

r = r(1+t0 : end - mod(length(r), M+Npx) + t0);

% Perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);
x_matrix = fft(r_matrix);

% Select useful samples
x_rcv = x_matrix(init_step:spacing:end, 1);
X_known = diag(ts)*sqrt(2);

% Compute G_est by dividing the received symbol by the transmitted one
G_est = x_rcv ./ ts;

% Solve LS for F*g=G_est where g is an 8x1 vector
F = dftmtx(M);
F = F(indices, 1:N2+1);
g_hat = (F' * F) \ (F' * G_est);
G_hat = fft(g_hat, M);


% Noise estimation
xhat = X_known * G_hat(init_step : spacing : end);
E = sum(abs(xhat - x_rcv).^2)/length(xhat);
est_sigma_w = E/M;


end
