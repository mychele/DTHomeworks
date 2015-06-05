function [G_hat, est_sigma_w] = OFDM_channel_estimation_2(snr, Npx, N2, t0)

%% Channel ESTIMATION for OFDM (second method)
% Send one block of data using 8 equally spaced groups of 4
% adjacent subchannels.

OFDM = true;
M = 512;
allowed_symb = 32;

block = ones(M, 1)*(-1-1i);
% Scale in order to double the power of tx symbols
ts = ts_generation(allowed_symb-1, 1) * sqrt(2);
sigma_ts = 4;


nsamples = 8;
symbpersegment = allowed_symb / nsamples;
indices = reshape(1:M, M/nsamples, nsamples);
indices = reshape(indices(1:symbpersegment, :), size(indices, 2)*symbpersegment, 1); % Second way
block(indices) = ts;

% Compute IDFT, add prefix, P/S
A = ifft(block);
A_pref = [A(end-Npx + 1:end); A];
s = reshape(A_pref, [], 1);

%% CHANNELIZATION

snr_lin = 10^(snr/10);
% Send over the noisy channel
[r, ~, ~] = channel_output(s, snr_lin, OFDM);


% --- Process at the receiver

r = r(1+t0 : end - mod(length(r), M+Npx) + t0);

% Perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);
x_matrix = fft(r_matrix);

% Select useful samples
x_rcv = x_matrix(indices, 1);

% Compute G_est by dividing the received symbol by the transmitted one
G_est = x_rcv ./ ts;

% Solve LS for F*g=G_est where g is an (N2+1)x1 vector
F = dftmtx(M);
F = F(indices, 1:N2+1);
g_hat = (F' * F) \ (F' * G_est);
G_hat = fft(g_hat, M);

% Noise estimation new method
est_sigma_w = 0;
for j=0:nsamples-1
   tempGest = G_est((j*symbpersegment + 1) : (j*symbpersegment + 1)+3);
   est_sigma_w = est_sigma_w + var(tempGest - mean(tempGest));
end
est_sigma_w = est_sigma_w / nsamples / M * sigma_ts;


end