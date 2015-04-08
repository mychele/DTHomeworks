close all
clear
clc

%% Load data
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
%z = randn(1000, 1);
K = length(z); % signal length


%% Spectral analysis

% AR model: order estimation
% compute variance of AR model and plot it to identify the knee
% it's computed up to K/5 - 1, don't know if it makes sense
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

% The knee is apparently at N = 3: use it as order of the AR model
plot_spectrum(z, 3);
ylim([-10 40])



%% Percentiles


D = round(44/0.228);
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