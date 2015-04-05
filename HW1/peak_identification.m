close all;
clear all;
clc;

% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
z = repmat(z, 40, 1);
K = length(z); % signal length

% TUNE THESE PARAMETERS AND SEE!
step = 64; %distance between the first two sample of each window
span = 192; %actual size of the window
overlap = span - step;

paddingfactor = 1;
padding = span * (paddingfactor - 1);  % Zero padding to apply to the windowed signal

locs_per = 1:K+padding; % just to initialize for the first iteration
peaks = zeros(K+padding,1); % just to initialize for the first iteration
acc_pks_per = zeros(span + padding, 1);

max_iter = floor((K-span)/step);

window_span = kaiser(span+padding, 16); %5.65);

h = zeros(span, 1);
h2 = h;

%figure
for i = 0:max_iter
    z_part = [z(i*step + 1: i*step + span); zeros(padding,1)];
    Z = fft(z_part.*window_span);
    periodogr = abs(Z).^2/span;
    [peaks, locs_per] = findpeaks(abs(periodogr));
    %plot(10*log10(periodogr))
    %pause
    edges = 0.5:paddingfactor:(span+padding+0.5);
    [N, ~] = histcounts(locs_per, edges);
    h = h + N.';
    %stem(locs_per/(span+padding), 10*log10(peaks)), ylim([-20, 40]), pause
    
    
    
    
    freqs = findSineNoise3(z_part, 32);
    %pause
    edges = (0.5:paddingfactor:(span+padding+0.5)) / (span+padding);
    [N, ~] = histcounts(freqs, edges);
    h2 = h2 + N.';
    
end

%normalize
h = h/max_iter;
h2 = h2/max_iter;

figure
bar((1:span)/span, h)
title('periodogram peaks accumulation')

figure
bar((1:span)/span, h2)
title('periodogram peaks accumulation 2')

disp(find(acc_pks_per > 0.7)/(span+padding))
