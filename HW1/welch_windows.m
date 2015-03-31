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

locs = 1:1000; % just to initialize for the first iteration
peaks = zeros(1000,1); % just to initialize for the first iteration


step = 13; %distance between the first two sample of each window
span = 130; %actual size of the window
overlap = span - step;
acc_locs = zeros(span, 1);

max_iter = floor((K-span)/step);

figure
for i = 0:max_iter
    z_part = z(i*step + 1: i*step + span);
    
    %% Compute different spectral analysis
    
    % PERIODOGRAM pg 84
    Z = fft(z_part);
    periodogr = abs(Z).^2/span;
    
    % compute WELCH estimator pg 85
    D = 50; % window size
    window = kaiser(D, 5.65);
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
    
    
    
    clear a autoc D fir_bs_1 N N_corr S upp_limit window window_correlogram
    
    %% Plot PSD estimate
    
    
%     plot((1:span)/span, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
%     hold on
%     plot((1:span)/span, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
%     hold on
%     plot((1:span)/span, 10*log10(periodogr), 'c:')
%     hold on
%     plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
%     hold on
%     plot(locs/span, 10*log10(peaks),'o','MarkerSize',12)
%     hold off
%     axis([0, 1, 0, 50])
%     legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)')
%     title(sprintf('Spectral analysis from sample %d', i*step))
    
    %pause();
    
    [peaks, locs] = findpeaks(P_welch);
    % 
    acc_locs(locs) = acc_locs(locs) + 1;
    
    
    
end

% normalize
acc_locs = acc_locs/i;

bar((1:span)/span, acc_locs)
axis([0, 1, 0, 0.4])
