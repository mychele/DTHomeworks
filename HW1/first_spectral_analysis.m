%% First spectral analysis
close all
clear
clc

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length

%% Spectral analysis

% AR model: order estimation.
% Compute variance of AR model and plot it to identify the knee.
% It's computed up to K/5 - 1
N_corr = ceil(K/5);
autoc = autocorrelation_biased(z, N_corr);
upp_limit = 30;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
    [~, sigma_w(N)] = arModel(N, autoc);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))
title('sigma_w of the AR model of the whole signal'), grid on
ylabel('sigma_w (dB)'), xlabel('Order of the AR(N) model')

% Plot the spectra of the signal according to various estimators
plot_spectrum(z, 3);
ylim([-10 40])

%% Peaks accumulation
% Find if a peak is present in more than one window of the signal

step = 64; % distance between the first two sample of each window
span = 96; % actual size of the window
overlap = span - step;
max_iter = floor((K-span)/step);

% Initialize
acc_locs_per = zeros(span, 1);
window_span = kaiser(span, 5.65);

for i = 0:max_iter
    z_part = z(i*step + 1: i*step + span);
    
    % PERIODOGRAM over span samples
    Z = fft(z_part.*window_span);
    periodogr = abs(Z).^2/span;
    
    % find local maxima
    [peaks, locs_per] = findpeaks(abs(periodogr));
    % accumulate local maxima
    acc_locs_per(locs_per) = acc_locs_per(locs_per) + 1;   
end

% normalize to the number of iterations
acc_locs_per = acc_locs_per/i;

figure
bar((1:span)/span, acc_locs_per)
axis([0, 1, 0, max(acc_locs_per)])
xlabel('Frequency')
title(sprintf('Histogram of periodogram peaks over %d windows', max_iter))

disp(find(acc_locs_per > 0.7)/span);

%% Whitening filter

% Pass the signal through a whitening filter in order to equalize its
% spectrum. The whitening filter will be obtained as the inverse of the AR 
% model filter. This is needed to identify peaks. As we do not care about 
% the phase in this stage, we can afford to pass it through something that
% is not with linear phase.

% Compute the AR model
[a, sigma_w] = arModel(3, autocorrelation(z, K/5));
[H, omega] = freqz(1, [1; a], K, 'whole');

% Plot the two frequency responses
figure
plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2));
hold on
plot(omega/(2*pi), 10*log10((1/sigma_w)*abs(1./H).^2));
hold off
legend('AR filter', 'Inverse of AR filter');
title('Frequency response of the AR filter and of its inverse');
axis([0 1 -40 40]);

% Plot the whitened result
figure
equalized_spectrum = abs((fft(z)).*(1./H)*(1/sigma_w));
plot((1:K-1)/K, 10*log10(equalized_spectrum(2:end))); hold on
title('Equalized spectrum of the given signal');
axis([0 1 0 20]);


%% Percentiles

D = round(200);
window = kaiser(D, 5.65);
S = round(D-K/64); %common samples
[Pm, PM, Psorted, ~] = findSine(z, window, S);

figure
plot((0:length(Pm)-1)/length(Pm), 10*log10(Pm).'), hold on
plot((0:length(PM)-1)/length(PM), 10*log10(PM).')
title('Minimum and maximum PSD across all windows')
ylabel('PSD (dB)'), ylim([-25 45])
% figure
% plot((0:length(PM)-1)/length(PM), 10*(log10(PM) - log10(Pm)).')
% title('Ratio between min and max PSD across all windows')

% Plot the percentiles
figure
hold on
percentileindices = round([1 30 70 99] * size(Psorted, 2) / 100);
percentileindices = max(percentileindices, ones(size(percentileindices)));
for i = percentileindices
    plot((0:length(Psorted)-1)/length(Psorted), 10*log10(Psorted(:, i)))
end
title('Percentiles of PSD (dB)')
legend('1', '30', '70', '99', 'Location', 'SouthEast')
ylim([-20 40])

%percentiles zoom
figure
hold on
percentileindices = round([1 10 20 30 70 80 90 99] * size(Psorted, 2) / 100);
percentileindices = max(percentileindices, ones(size(percentileindices)));
for i = percentileindices
    plot((0:length(Psorted)-1)/length(Psorted), 10*log10(Psorted(:, i)))
end
title('Percentiles of PSD (dB)')
legend('1', '10', '20', '30', '70', '80', '90', '99', 'Location', 'SouthEast')
ylim([-20 20]), xlim([.6 .9])