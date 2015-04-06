close all;
clear all;
clc;

%% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length



step = 20; %distance between the first two sample of each window
span = 100; %actual size of the window
overlap = span - step;

locs_w = 1:1000;
locs_c = 1:1000;
locs_per = 1:1000; % just to initialize for the first iteration
peaks = zeros(1000,1); % just to initialize for the first iteration
peaks_w = zeros(1000,1);
peaks_c = zeros(1000,1);
acc_locs_per = zeros(span, 1);
acc_locs_w = zeros(span, 1);
acc_locs_c = zeros(span, 1);

max_iter = floor((K-span)/step);

window_span = kaiser(span, 5.65);

figure
for i = 0:max_iter
    z_part = z(i*step + 1: i*step + span);
    
    %% Compute different spectral analysis
    
    % PERIODOGRAM pg 84
    Z = fft(z_part.*window_span);
    periodogr = abs(Z).^2/span;
    
    % compute WELCH estimator pg 85
    D = 50; % window size
    window = kaiser(D, 16);
    S = D/2; %common samples
    P_welch = welchPsd(z_part, window, S);
    
    % CORRELOGRAM
    N_corr = ceil(span/5); % N_corr is the order of the autocorrelation estimate
    window_correlogram = kaiser(2*N_corr + 1, 5.65); % window centered around N_corr
    correlogram = correlogramPsd(z_part, window_correlogram, N_corr);
    
    
    % AR(3)
    autoc = autocorrelation(z_part, N_corr);
    N = 3;
    [a, sigma_w] = arModel(N, autoc);
    [H, omega] = freqz(1, [1; a], span, 'whole');
    
    
    
    clear a autoc D fir_bs_1 N N_corr S upp_limit window %window_correlogram
    
    %% Plot PSD estimate
    
    
%     plot((1:span)/span, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
%     hold on
%     plot((1:span)/span, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
%     hold on
%     plot((1:span)/span, 10*log10(periodogr), 'c:')
%     hold on
%     plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
%     hold on
%     plot(locs_w/span, 10*log10(peaks_w),'o','MarkerSize',12)
%     hold off
%     axis([0, 1, 0, 50])
%     legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)')
%     title(sprintf('Spectral analysis from sample %d', i*step))
%     
%     pause();
    
    [peaks, locs_per] = findpeaks(abs(periodogr));
    [peaks_w, locs_w] = findpeaks(P_welch);
    [peaks_c, locs_c] = findpeaks(abs(correlogram));

    % 
    acc_locs_per(locs_per) = acc_locs_per(locs_per) + 1;
    acc_locs_w(locs_w) = acc_locs_w(locs_w) + 1;
    acc_locs_c(locs_c) = acc_locs_c(locs_c) + 1;
    
    
end

% normalize
acc_locs_per = acc_locs_per/i;
acc_locs_w = acc_locs_w/i;
acc_locs_c = acc_locs_c/i;


figure
% subplot(1, 3, 1)
bar((1:span)/span, acc_locs_per)
axis([0, 1, 0, max(acc_locs_per)])
title('periodogram peaks accumulation')
% subplot(1, 3, 2)
% bar((1:span)/span, acc_locs_w)
% axis([0, 1, 0, max(acc_locs_per)])
% title('welch')
% subplot(1, 3, 3)
% bar((1:span)/span, acc_locs_c)
% axis([0, 1, 0, max(acc_locs_per)])
% title('correlogram')

disp(find(acc_locs_per > 0.7)/span)
