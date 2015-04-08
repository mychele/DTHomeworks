close all;
clear all;
clc;

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
%z = randn(1000, 1);
K = length(z); % signal length


%% Compute different spectral analysis

% PERIODOGRAM pg 84
Z = fft(z);
periodogr = abs(Z).^2/K;

% compute WELCH estimator pg 85 
D = 200; % window size
window = kaiser(D, 5.65);
S = D/2; %common samples
P_welch = welchPsd(z, window, S);

% CORRELOGRAM
N_corr = ceil(K/5); % N_corr is the order of the autocorrelation estimate
window_correlogram = kaiser(2*N_corr + 1, 5.65); % window centered around N_corr
correlogram = correlogramPsd(z, window_correlogram, N_corr);

% AR model: order estimation
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
autoc = autocorrelation(z, N_corr);
upp_limit = 60;
sigma_w = zeros(1, upp_limit);
for N = 1:upp_limit
    [~, sigma_w(N)] = arModel(N, autoc);
end
figure, plot(1:upp_limit, 10*log10(sigma_w))

% the knee is apparently at N = 3
% compute the vector of coefficients a
N = 3;
[a, sigma_w] = arModel(N, autoc);
[H, omega] = freqz(1, [1; a], K, 'whole');

clear a autoc D fir_bs_1 N N_corr S upp_limit window window_correlogram



%% Plot PSD estimate

figure
plot((1:K)/K, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
hold on
plot((1:K)/K, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
hold on
plot((1:K)/K, 10*log10(periodogr), 'c:')
hold on
plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
hold off
axis([0, 1, -10, 40])
legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)', 'Location', 'SouthEast')
title('Spectral analysis')



%% Percentiles


D = round(44/0.228);
window = kaiser(D, 5.65);
S = round(D-K/64); %common samples
[Pm, PM, Psorted, ~] = findSine(z, window, S);

%figure, plot(10*log10(abs(fft_mean)));

figure
plot((0:length(Pm)-1)/length(Pm), 10*log10(Pm).'), hold on
plot((0:length(PM)-1)/length(PM), 10*log10(PM).')
title('Minimum and maximum PSD across all windows')
ylabel('PSD (dB)')
figure
plot((0:length(PM)-1)/length(PM), 10*(log10(PM) - log10(Pm)).')
title('Ratio between min and max PSD across all windows')

% percentiles
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




% fft_meanmean = zeros(size(z));
% Dvalues = 100:4:200;
% for D=Dvalues
%     window = kaiser(D, 5.65);
%     S = round(D*0.9); %common samples
%     [~, ~, ~, fft_mean] = findSine(z, window, S);
%     fft_meanmean = fft_meanmean + fft_mean;
% end
% figure, plot(10*log10(abs(fft_meanmean / length(Dvalues))));