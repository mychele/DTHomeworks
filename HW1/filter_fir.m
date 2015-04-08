close all
clear
clc

%% Load data
load('data for hw1.mat');
z = z.'; % make a column vector
load('hp18.mat');
hp18 = hp18.'; % make a column vector
K = length(z); % signal length

%% Complex BPF
% --- Compute the coefficients
f0 = 0.770; % estimated by inspection on the PSD + confirmed by method in ampphase_estimation_rls
% cfirpm has a strange behaviour, the center of the band is -(1-f0)*2
freq_delimiters = [f0 - 0.02, f0 - 0.002, f0 + 0.002, f0 + 0.02]; % limit of don't care regions, left and right of f0
matlab_correct_setting = -2*(1-freq_delimiters);
% bandpass filter designed with cfirmpm
bpf = cfirpm(58, [-1, matlab_correct_setting, 1], @bandpass);

% --- Plot frequency response of bandpass filter
% DTFTplot(bpf, 50000);
% ylim([-40 0])
% title('Freq resp of the BPF filter')

%% Filter the signal with HPF + BPF
linesfilter = conv(hp18, bpf);
% normalize linesfilter
linesfilter = linesfilter / max(abs(fft([linesfilter zeros(1, 5000 - length(linesfilter))])));
N_linesfilter = length(linesfilter) - 1; % Order of the filter
z_lines = filter(linesfilter, 1, z);
DTFTplot(linesfilter, 10000); % Plot filter's freq resp
title('Freq response of HPF + BPF (dB)')
ylim([-50 0]), grid on
% Compensate delay (still half transient left)
z_lines = z_lines( (N_linesfilter/2 + 1) : length(z_lines));


%% Compute complementary filter and get "continuous PSD" part

% Compute the complementary of the filter we just used
linesfilter_compl = -linesfilter;
linesfilter_compl(N_linesfilter/2 + 1) = linesfilter_compl(N_linesfilter/2 + 1) + 1;
DTFTplot(linesfilter_compl, 10000);
title('Freq response of complementary filter (dB)')
ylim([-25 5]), grid on

% Filter original signal
z_continuous = filter(linesfilter_compl, 1, z);
% Compensate delay (still half transient left)
z_continuous = z_continuous( (N_linesfilter/2 + 1) : length(z_continuous));


%% Export the two signals

save('split_signal', 'z_continuous', 'z_lines');

%% See that diff is zero

delayonly = [zeros(N_linesfilter/2, 1); 1];
delayed_z = filter(delayonly, 1, z);
% Compensate delay
delayed_z = delayed_z( (N_linesfilter/2 + 1) : length(delayed_z));
diff = delayed_z - (z_continuous + z_lines);
disp(['Max magnitude of the difference between the original signal and the ', ...
    'sum of its two components (lines and continuous) is ', num2str(max(abs(diff)))])

%% Remove mean and plot spectral analysis (without AR models)

% Remove mean from continuous and lines part
z_continuous = z_continuous - mean(z_continuous);
z_lines = z_lines - mean(z_lines);

% Plot spectral analysis of the signal and its two components (without AR models)
plot_spectrum(z_lines, 0); 
axis([0 1 -10 40]), title('Spectral analysis of spectral line part')
plot_spectrum(z_continuous, 0);
axis([0 1 -10 40]), title('Spectral analysis of continuous part')
plot_spectrum(z, 3);
axis([0 1 -10 40]), title('Spectral analysis of original signal')



%% AR model for continuous part

% Find the knee of sigma_w
N_corr = floor(length(z_continuous)/5);
autoc_cont = autocorrelation_biased(z_continuous, N_corr);
upp_limit = 30;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
   [~, sigma_w(N)] = arModel(N, autoc_cont);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))
title('sigma_w of the AR model of the continuous part'), grid on
ylabel('sigma_w (dB)'), xlabel('Order of the AR(N) model')

% Choose order for AR and compute the vector of coefficients a
N = 3;
[a_cont, sigma_w_cont] = arModel(N, autoc_cont);
[H, omega] = freqz(1, [1; a_cont], K, 'whole');
figure, plot(omega/(2*pi), 10*log10(sigma_w_cont*abs(H).^2), 'Color', 'm', 'LineWidth', 1)
axis([0, 1, -10, 40])
title('AR model of the continuous part')
xlabel('Normalized frequency'), ylabel('Magnitude (dB)')

figure, zplane(roots([1;a_cont]))
%% AR model for spectral lines

% Find the knee of sigma_w
N_corr = floor(length(z_lines)/5);
autoc = autocorrelation_biased(z_lines, N_corr);
upp_limit = 30;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
   [~, sigma_w(N)] = arModel(N, autoc);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))
title('sigma_w of the AR model of the spectral lines'), grid on
ylabel('sigma_w (dB)'), xlabel('Order of the AR(N) model')

% Choose order for AR and compute the vector of coefficients a
N = 5;
[a_lines, sigma_w_lines] = arModel(N, autoc);
[H, omega] = freqz(1, [1; a_lines], K, 'whole');
figure, plot(omega/(2*pi), 10*log10(sigma_w_lines*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
axis([0, 1, -35, 15])
title('AR model of the spectral lines')
xlabel('Normalized frequency'), ylabel('Magnitude (dB)')

figure, zplane(roots([1;a_lines]))
