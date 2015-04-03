close all;
clear all;
clc;

% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length

step = 50; %distance between the first two sample of each window
span = 200; %actual size of the window
overlap = span - step;

padding = 0;  % Zero padding to apply to the windowed signal

locs_per = 1:K+padding; % just to initialize for the first iteration
peaks = zeros(K+padding,1); % just to initialize for the first iteration
acc_locs_per = zeros(span + padding, 1);

max_iter = floor((K-span)/step);

window_span = kaiser(span+padding, 5.65);

frequencies = zeros(0,0);

figure
for i = 0:max_iter
    z_part = [z(i*step + 1: i*step + span); zeros(padding,1)];
    
    Z = fft(z_part.*window_span);
    periodogr = abs(Z).^2/span;
    
    plot(10*log10(periodogr))
    pause(0.5)
        
    [peaks, locs_per] = findpeaks(abs(periodogr));
    
    frequencies = [frequencies, locs_per.'];

    acc_locs_per(locs_per) = acc_locs_per(locs_per) + 1;
    
end

normalize
acc_locs_per = acc_locs_per/i;

figure
h = histogram(frequencies/(span+padding), span+padding);
bar((1:(span+padding))/(span+padding), acc_locs_per)
axis([0, 1, 0, max(acc_locs_per)])
title('periodogram peaks accumulation')

disp(find(acc_locs_per > 0.7)/(span+padding))
