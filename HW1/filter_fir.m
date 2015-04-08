close all
clear
clc

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
hp18 = load('hp18.mat');
hp18 = hp18.hp18.'; % make a column vector
%z = z - mean(z); % remove average
K = length(z); % signal length
autoc = autocorrelation(z, length(z)/5);


%% BPF plus complex bandpass filter
% --- Compute the coefficients
f0 = 0.770;
% cfirpm has a strange behaviour, the center of the band is -(1-f0)*2
freq_delimiters = [f0 - 0.02, f0 - 0.002, f0 + 0.002, f0 + 0.02]; % limit of don't care regions, left and right of f0
matlab_correct_setting = -2*(1-freq_delimiters);
bpf = cfirpm(58, [-1, matlab_correct_setting, 1], @bandpass);

% --- Plot frequency response of bandpass filter
DTFTplot(bpf, 50000);
ylim([-40 0])
title('Freq resp of the BPF filter')

%% Filter the signal with HPF + BPF
linesfilter = conv(hp18, bpf);
N_linesfilter = length(linesfilter) - 1; % Order of the filter
z_lines = filter(linesfilter, 1, z);
DTFTplot(linesfilter, 10000); % Plot filter's freq resp
title('Freq response of HPF + BPF (dB)')
% Discard transient
z_lines = z_lines( (N_linesfilter/2 + 1) : length(z_lines));


%% Compute complementary filter and get "continuous PSD" part

% Compute the complementary of the filter we just used
linesfilter_compl = -linesfilter;
linesfilter_compl(N_linesfilter/2 + 1) = linesfilter_compl(N_linesfilter/2 + 1) + 1;
DTFTplot(linesfilter_compl, 10000);
title('Freq response of complementary filter (dB)')

% Filter original signal
z_continuous = filter(linesfilter_compl, 1, z);
% Discard transient
z_continuous = z_continuous( (N_linesfilter/2 + 1) : length(z_continuous));


%% Export the two signals

save('split_signal', 'z_continuous', 'z_lines');



%% Plots and stuff

% See that diff is zero
delayonly = [zeros(N_linesfilter/2, 1); 1];
delayed_z = filter(delayonly, 1, z);
% Discard transient
delayed_z = delayed_z( (N_linesfilter/2 + 1) : length(delayed_z));
diff = delayed_z - (z_continuous + z_lines);
disp(['Max magnitude of the difference between the original signal and the ', ...
    'sum of its two components (lines and continuous) is ', num2str(max(abs(diff)))])

z_continuous = z_continuous - mean(z_continuous);
% Plot original signal and continuous and line (without AR models)
plot_spectrum(z_lines, 0); 
ylim([-10 40]), title('Spectral analysis of spectral line part')
plot_spectrum(z_continuous, 0);
ylim([-10 40]), title('Spectral analysis of continuous part')
plot_spectrum(z, 4);
ylim([-10 40]), title('Spectral analysis of original signal')
%% AR model for continuous part


% --- Recompute everything for AR
% remove deterministic components

% Find the knee of sigma_w
N_corr = floor(length(z_continuous)/5);
autoc_cont = autocorrelation(z_continuous, N_corr);
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
   [~, sigma_w(N)] = arModel(N, autoc_cont);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))  % Uncomment to plot sigma

% Choose order for AR and compute the vector of coefficients a
N = 4;
[a, sigma_w] = arModel(N, autoc_cont);
[H, omega] = freqz(1, [1; a], K, 'whole');
figure, plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
axis([0, 1, -10, 40])

%% AR model for spectral lines

% --- Recompute everything for AR

% Find the knee of sigma_w
N_corr = floor(length(z_lines)/5);
autoc = autocorrelation(z_lines, N_corr);
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
   [~, sigma_w(N)] = arModel(N, autoc);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))  % Uncomment to plot sigma

% Choose order for AR and compute the vector of coefficients a
N = 2;
[a, sigma_w] = arModel(N, autoc);
[H, omega] = freqz(1, [1; a], K, 'whole');
figure, plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
axis([0, 1, -40, 10])



