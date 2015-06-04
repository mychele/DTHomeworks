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

ts = ts_generation(allowed_symb-1, 1) * sqrt(2);
init_step = 1; % < 16
indices = init_step : spacing : init_step + spacing*(allowed_symb-1);
block(indices) = ts;

% Compute IDFT, add prefix, P/S
A = ifft(block);
A_pref = [A(end-Npx + 1:end); A];
s = reshape(A_pref, [], 1);

%% CHANNELIZATION

snr = 6; %dB
snr_lin = 10^(snr/10);
% Send over the noisy channel
[r, sigma_w, g] = channel_output(s, snr_lin, OFDM);
g = g(1+t0 : end);   % Take t0 into account (just to plot stuff)
G = fft(g, 512);
G = G(:);


%% Process at the receiver

r = r(1+t0 : end - mod(length(r), M+Npx) + t0);

% Perform the DFT
r_matrix = reshape(r, M+Npx, []);
r_matrix = r_matrix(Npx + 1:end, :);
x_matrix = fft(r_matrix);

% Select useful samples
x_rcv = x_matrix(indices, 1);
x_known = diag(ts);

% Compute G_est by dividing the received symbol by the transmitted one
G_est = x_rcv ./ ts;

% Solve LS for F*g=G_est where g is an 8x1 vector: do it for different values of N2
F_complete = dftmtx(M);

for N2 = 1:Npx
   F = F_complete(indices, 1:N2+1);
   g_hat = (F' * F) \ (F' * G_est);
   g_est = ifft(G_est);
   G_hat = fft(g_hat, M);
   
   % Noise estimation
   xhat = x_known * G_hat(indices);
   E = sum(abs(xhat - x_rcv).^2)/length(xhat);
   est_sigma_w(N2) = E/M; %#ok<SAGROW>
end



%% Plot

figure, hold on
plot(1:Npx, 10*log10(est_sigma_w))
plot([1 Npx], 10*log10(sigma_w)*[1 1])
title('Error functional for channel estimation varying N2')
xlabel('N2'), ylabel('Error functional (dB)')
grid on, box on