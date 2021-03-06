close all;
clear all;
clc;

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
hp18 = load('hp18.mat');
hp18 = hp18.hp18.'; % make a column vector
%z = z - mean(z); % remove average
K = length(z); % signal length
autoc = autocorrelation(z, length(z)/5);


%% HPF plus Wiener filter pg 137 BC
N_wien = 31;        % Number of coefficients of the Wiener filter
w0 = 2*pi*0.78;     % 0.78 is probably our freq
E_w0 = exp(1i * w0 * (0 : N_wien-1)).';
A = sqrt(30); % see plot of real autoc, rule of thumb (actually doesn't change much)
B=A;
snr = A^2/autoc(1); % autoc(1) should be sigma^2 of the noise, if white
gain = B/A;
D = 100;    % delay

% --- Compute the coefficients
c_opt = gain * snr * exp(- 1i * w0 * D) * E_w0 / (1 + N_wien*snr);


% --- Plot frequency response of Wiener filter
DTFTplot(c_opt, 50000);
ylim([-20 0])
title('Freq resp of Wiener filter')



%% Filter the signal with LPF + Wiener
linesfilter = conv(hp18, c_opt);
N_linesfilter = length(linesfilter) - 1; % Order of the filter
z_lines = filter(linesfilter, 1, z);
DTFTplot(linesfilter, 10000); % Plot filter's freq resp
title('Freq response of LPF + Wiener filter (dB)')
% Discard transient
z_lines = z_lines( (N_linesfilter/2 + 1) : length(z_lines));



%% Analysis

% --- Plot spectrum of the spectral lines
plot_spectrum(z_lines, 1); % 1 is the order of the desired AR model
title('Spectral analysis of the signal after filtering'), legend('Location', 'NorthWest')


% --- Recompute everything for AR

% Find the knee of sigma_w
%N_corr = length(spectral_lines)/5;
%autoc = autocorrelation(spectral_lines, N_corr);
%upp_limit = 60;
%sigma_w = zeros(1, upp_limit);
%for N = 1:upp_limit
%    [~, sigma_w(N)] = arModel(N, autoc);
%end
%figure, plot(1:upp_limit, 10*log10(sigma_w))  % Uncomment to plot sigma

% Choose order for AR and compute the vector of coefficients a
%N = 1;
%[a, sigma_w] = arModel(N, autoc);
%[H, omega] = freqz(1, [1; a], K, 'whole');
%figure, plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
%axis([0, 1, -40, 10])




%% Compute complementary filter and get "continuous PSD" part

% Compute the complementary of the filter we just used
linesfilter_compl = -linesfilter;
linesfilter_compl(N_linesfilter/2) = linesfilter_compl(N_linesfilter/2) + 1;
DTFTplot(linesfilter_compl, 10000);
title('Freq response of complementary filter (dB)')

% Filter original signal
z_continuous = filter(linesfilter_compl, 1, z);
% Discard transient
z_continuous = z_continuous( (N_linesfilter/2 + 1) : length(z_continuous));


%% Plots and stuff

% Plot original signal and continuous (without AR models)
plot_spectrum(z_continuous, 0);
ylim([-10 40]), title('Spectral analysis of continuous part')
plot_spectrum(z, 0);
ylim([-10 40]), title('Spectral analysis of original signal')

% See that diff is zero
delayonly = [zeros(N_linesfilter/2 - 1, 1); 1];
delayed_z = filter(delayonly, 1, z);
% Discard transient
delayed_z = delayed_z( (N_linesfilter/2 + 1) : length(delayed_z));
diff = delayed_z - (z_continuous + z_lines);
disp(['Max magnitude of the difference between the original signal and the ', ...
    'sum of its two components (lines and continuous) is ', num2str(max(abs(diff)))])

%% Export the two signals

save('split_signal', 'z_continuous', 'z_lines');