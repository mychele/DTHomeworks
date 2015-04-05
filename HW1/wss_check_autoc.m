close all
clear all
clc

%% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
K = length(z); % signal length

step = 100; %distance between the first two sample of each window
span = 100; %actual size of the window
N_corr = ceil(span/5);

figure
max_iter = floor((K-span)/step);
for i = 0:max_iter
    z_part = z(i*step + 1: i*step + span);
    
    autoc = autocorrelation(z_part, N_corr);
    
    plot(1:N_corr + 1, real(autoc), 'r')
    hold on
    plot(1:N_corr + 1, imag(autoc), 'c')
    
    
    %pause();
    
end
hold off
xlabel('lag')
ylabel('autocorrelation')
legend('real(autocorr)', 'imag(autocorr)')