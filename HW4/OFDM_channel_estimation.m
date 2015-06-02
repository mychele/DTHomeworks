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
block = ones(M, 1)*(-1-1i);

%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%
% NOT SURE IF IT'S OK TO SCALE IN THIS WAY
% Remember: for the symbols on which the estimation is performed, for a
% given snr (computed with the usual sigma_a^2 = 2), the power of the
% "estimation symbols" is doubled (-> better snr)
% Note that the variance of the noise at the receiver, after the DFT, is
% multplied by M, therefore it could be high
% Scale in order to double the power of tx symbols
ts = ts_generation(31, 1);
init_step = 1; % < 16
block(init_step:spacing:end) = ts * sqrt(2);

A = ifft(block);
A_pref = [A(end-Npx + 1:end); A];

s = reshape(A_pref, [], 1);

%% CHANNELIZATION
snr = 6; %dB
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
x_rcv = x_matrix(init_step:spacing:end, 1)/sqrt(2);
x_rcv(length(x_rcv) + 1) = x_rcv(1);
X_known = diag(ts_generation(31, 2));

phi = X_known'*X_known;
theta = X_known'*x_rcv;

G_est = phi \ theta;

% InterpolaTION
f = init_step:spacing:M+init_step;
f_fine = 1:1:M+init_step*spacing;
G_est_plusone = interp1(f, G_est, f_fine, 'spline');

% Drop the samples outside one period
x_rcv = x_rcv(1:end-1);
X_known = diag(ts_generation(31, 1));
G_est = G_est(1:end-1);
G_est_complete = G_est_plusone(1:end-init_step*spacing);
f_support = 1:length(G_est_complete);

figure, 
subplot 211
plot(real(G_est_complete)), hold on
plot(init_step:spacing:M, real(G_est), '-h'), hold on
plot(real(G)),
title(strcat('Comparison between estimated - LS+interpol - and real at ', num2str(snr), ' dB'))
legend('real(G_{est}) iterpolated', 'real(G_{est}) not interpolated', 'real(G)'), xlabel('i - subchannels'), ylabel('Real(G)'),
grid on
subplot 212
plot(imag(G_est_complete)), hold on
plot(init_step:spacing:M, imag(G_est), '-h'), hold on
plot(imag(G)),
title(strcat('Comparison between estimated - LS+interpol - and real at ', num2str(snr), ' dB'))
legend('imag(G_{est}) iterpolated', 'imag(G_{est}) not interpolated', 'imag(G)'), xlabel('i - subchannels'), ylabel('imag(G)'),
grid on


% Estimation of the noise
xhat = X_known*sqrt(2)*G_est;
x_rcv = x_matrix(1:spacing:end);
E = sum(abs(xhat - x_rcv).^2)/length(xhat);
est_sigma_w = E/M;
fprintf('Est sigma_w^2 = %d\n', est_sigma_w);
fprintf('Real sigma_w^2 = %d', sigma_w);

