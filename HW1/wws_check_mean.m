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
span = step; %actual size of the window
N_corr = ceil(span/5);

figure
max_iter = floor((K-span)/step);

step_short = step/5;
span_short = span/5;

max_iter_in = floor((span - span_short)/step_short);
for i = 0:max_iter
    z_part = z(i*step + 1: i*step + span);
    
    z_part_mean(i+1) = mean(z_part); % remove average
    
    
    for r = 0:max_iter_in
        % note that since the samples are really few the variance of the
        % estimator is high
        % it doesn't make sense compute the estimate of the
        % autocorrelation, since there are too few samples
        z_in = z_part(r*step_short + 1: r*step_short + span_short);
        z_part_in(i+1, r+1) = mean(z_in);
    end
    
end

% comment to study autocorrelation
figure(1)
plot(real(z_part_mean))
hold on
plot(imag(z_part_mean), 'r')
hold on
plot(max(real(z))*ones(max_iter + 1, 1), 'c')
hold on
plot(min(real(z))*ones(max_iter + 1, 1), 'c')
hold on
plot(max(imag(z))*ones(max_iter + 1, 1), 'y')
hold on
plot(min(imag(z))*ones(max_iter + 1, 1), 'y')
xlabel('window')
ylabel('mean and max min')
legend('real(mean)', 'imag(mean)', 'max(real(z))', 'min(real(z))', 'max(imag(z))', 'min(imag(z))')


for i = 1:max_iter+1
    figure
    plot(real(z_part_in(i, :)))
    hold on
    plot(imag(z_part_in(i, :)), 'r')
    hold on
    plot(max(real(z))*ones(max_iter_in + 1, 1), 'c')
    hold on
    plot(min(real(z))*ones(max_iter_in + 1, 1), 'c')
    hold on
    plot(max(imag(z))*ones(max_iter_in + 1, 1), 'y')
    hold on
    plot(min(imag(z))*ones(max_iter_in + 1, 1), 'y')
    xlabel('window')
    ylabel('mean and max min')
    legend('real(mean)', 'imag(mean)', 'max(real(z))', 'min(real(z))', 'max(imag(z))', 'min(imag(z))')
end