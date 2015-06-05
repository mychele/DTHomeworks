%% Channel ESTIMATION for OFDM (second method)
% Send one block of data using 8 equally spaced groups of 4
% adjacent subchannels.

%#ok<*SAGROW>

clear, close all
numsim = 500;
OFDM = true;
M = 512;
allowed_symb = 32;
Npx = 7;
N2 = 4;
t0 = 5;

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


snr = 6; %dB
snr_lin = 10^(snr/10);
% Send over the noisy channel
[r, sigma_w, g] = channel_output(s, snr_lin, OFDM);
g = g(1+t0 : end);   % Take t0 into account (just to plot stuff)
G = fft(g, 512);
G = G(:);


% --- Process at the receiver

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

% Solve LS for F*g=G_est where g is an 8x1 vector
F = dftmtx(M);
F = F(indices, 1:N2+1);
g_hat = (F' * F) \ (F' * G_est);
g_est = ifft(G_est);
G_hat = fft(g_hat, M);

% Noise estimation new method
est_sigma_w = 0;
for j=0:nsamples-1
   tempGest = G_est((j*symbpersegment + 1) : (j*symbpersegment + 1)+3);
   est_sigma_w = est_sigma_w + var(tempGest - mean(tempGest));
end
est_sigma_w = est_sigma_w / nsamples / M * sigma_ts;


% Error on the estimate of G
est_err = sum(abs(G_hat - G).^2);



%% Plots

fprintf('Est sigma_w^2 = %.2d\n', est_sigma_w);
fprintf('Real sigma_w^2 = %.2d\n', sigma_w);

figure, hold on
stem(0:Npx, abs(g))
stem(0:N2, abs(g_hat), 'x')
stem(0:10, abs(g_est(1:11)), '^')
legend('Actual g', 'g_{hat}', 'IDFT of G_{est}')
grid on, box on

figure,
subplot 211
plot(real(G)), hold on
plot(real(G_hat))
plot(indices, real(G_est), '^')
title(strcat('Real part of true and estimated G at ', num2str(snr), ' dB'))
legend('real(G)', 'real(G_{hat})', 'real(G_{est})'), xlabel('i - subchannels'), ylabel('Real(G)'),
grid on, xlim([1, M])

subplot 212
plot(imag(G)), hold on
plot(imag(G_hat))
plot(indices, imag(G_est), '^')
title(strcat('Real part of true and estimated G at ', num2str(snr), ' dB'))
legend('imag(G)', 'imag(G_{hat})', 'imag(G_{est})'), xlabel('i - subchannels'), ylabel('imag(G)'),
grid on, xlim([1, M])