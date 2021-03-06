%% Clean up and initialize useful quantities
clear
close all
clc
rng default

%% Data initialization

% The data_init script initializes all the input parameters that don't 
% change across the simulation
data_init;

% Plot residual energy of IIR filter in order to determine transient length
[h_dopp, ~] = impz(b_dopp, a_dopp);
residual_nrg = sum(h_dopp.^2) - cumsum(h_dopp.^2);
plot(0:length(residual_nrg)-1, 10*log10(residual_nrg)), 
xlim([0 100]), grid on, box on, title('Residual energy of I.R. of the IIR filter')
xlabel('Transient length (samples)'), ylabel('Residual energy (dB)')
% Doppler filter shape
[Hf, f] = freqz(b_dopp, a_dopp, 1000, 'whole');
figure, plot(f/(2*pi), 20*log10(abs(Hf)))
ylim([0, 12])
title('Frequency response of the Doppler filter')

%% Criterion for N_h
% We consider the error that is introduced by dropping the tail of the of
% the sampled exp function that describes the aleatory part of the PDP. The
% entire function is normalized to sum to 1-C^2.

M_complete = 1/tau_rms*exp(-(0:894)*Tc/tau_rms);
M_complete = M_complete.*(1-C^2)/sum(M_complete);

for N_h = 1:10
    % Only consider the tail
    deltaM = M_complete(N_h+1:end);
    lambda_n(N_h) = 1/(snr_lin * sum(abs(deltaM)));
end

figure
plot(1:length(lambda_n), 10*log10(lambda_n), 'd-')
grid on, title('\Lambda_n')
xlabel('N_h'), ylabel('\Lambda_n [dB]')
axis([1 4 -2 10]), ax = gca; ax.XTick = 1:5;

%% Display the PDP for the channel

% The PDP is the sampling of a continuous time exponential PDP with
% tau_rms / T = 0.3, so tau_rms / Tc = 1.2
% (See 4.224 for reference)
N_h = 3; % As determined before
tau = 0:Tc:N_h-1;
M_iTc = 1/tau_rms * exp(-tau/tau_rms); % 
% normalize pdp: it must be sum(E[|htilde_i|^2]) = 1 - C^2
M_iTc = M_iTc.*(1-C^2)/sum(M_iTc);
M_d = sum(M_iTc);

pdp = M_iTc;

% add LOS component power
pdp(1) = pdp(1) + C^2;
pdp_log = 10*log10(pdp); 

% Plot the normalised PDP
figure, stem(tau, pdp_log), title('PDP'), xlabel('iTc'), ylabel('E[|h_i(nTc)|^2]');
grid on;
axis([-0.25 2.25 -15 0]), ax = gca; ax.XTick = 0:2;
legend('PDP', 'Location', 'SouthWest')

%% Generation of impulse responses
% The code to generate the impulse responses is externalized in
% channel_generator script, which will be invoked both in the first and
% second exercise
channel_generator;

%% Show the behavior of |h_1(nTc)| for n = 0:1999, dropping the transient

figure, hold on
plot(0:1999, abs(h_mat(2, 1:2000).'))
grid on, box on, xlabel('nT_c'), ylabel('|h_1(nT_c)|')
title('|h_1(nT_C)|')

%% Show all |h_i|'s
% Just for debugging purposes

figure, hold on
plot(0:9999, abs(h_mat(:, 1:10000).'))
grid on, box on, xlabel('Time samples'), ylabel('|h_i(nT_C)|_{dB}')
legend('h_0', 'h_1', 'h_2')
title('|h_i|')

%% Histogram of h_1

% Plot of the required histogram
figure, histogram(abs(h_mat(2, 1:1000)).'/sqrt(M_iTc(2)), 20, ...
    'Normalization','pdf', 'DisplayStyle', 'stairs');
title('Experimental PDF from 1000 samples of hbar_1')
xlabel('hbar_1');

% This plot can be used to explain that because of the correlation we can't
% get a nice pdf (too little samples, correlation in peaks)
figure, 
subplot 121
histogram(abs(h_mat(2, 1:1000).')/sqrt(M_iTc(2)), 20, ...
    'Normalization','pdf', 'DisplayStyle', 'stairs');
camroll(90)
title('1000 samples of hbar_1')
subplot 122
plot(0:999, abs(h_mat(2, 1:1000).')/sqrt(M_iTc(2)));
title('Realization of hbar_1 over which the histogram is computed');
xlabel('Samples');
ylabel('hbar_1');
grid on

% Here we show that the experimental PDF gets better with more samples
figure
histogram(abs(h_mat(2, 1: 100000).')/sqrt(M_iTc(2)), 20, ...
    'Normalization','pdf', 'DisplayStyle', 'stairs')
title('100000 samples of hbar_1 vs Rayleigh pdf')
hold on
a = 0:0.01:3;
plot(a, 2.*a.*exp(-a.^2), 'LineWidth', 1.5);  % Theoretical PDF (page 308, BC)
hold off
legend('hbar_1', 'Rayleigh pdf');
ylabel('p_{hbar_1(kT_C)}(a)')
xlabel('a');

%% Simulation in order to compute the histogram of |h1(151Tc)|/sqrt(E(|h1(151Tc)|^2))
% This simulation repeats for numexp times, indipendently, the generation
% of the impulse response for ray 1. The task is the same as in the more
% general channel_generator. However, in order to speed up the simulation
% and since only one ray is involved, we don't invoke channel_generator
% script as before but generate the required impulse response only.
% Moreover, since we are interested in the 151th sample after the transient
% we can generate shorter impulse responses at each iteration.
numsim = 1000;
h_samples_needed = 200000 + ceil(Tp/Tc*length(h_dopp)); 
% Some will be dropped because of transient, since
% enough time, memory and computational power are available 
w_samples_needed = ceil(h_samples_needed / Tp);
transient = ceil(Tp/Tc*length(h_dopp));

h_1 = zeros(numsim, 1);
for k = 1:numsim
    disp(k)
    w = wgn(w_samples_needed,1,0,'complex');
    hprime = filter(b_dopp, a_dopp, w);
    % Interpolation
    t = 1:length(hprime);
    t_fine = Tq/Tp:Tq/Tp:length(hprime);
    h_fine = interp1(t, hprime, t_fine, 'spline');
    % Drop the transient and energy scaling
    h_notrans = h_fine(50000:end)*sqrt(M_iTc(2));
    % Energy scaling
    h_1(k) = h_notrans(152);  
end

figure, 
histogram(abs(h_1)/sqrt(sum(abs(h_1).^2)/length(h_1)), 20, ...
    'Normalization','pdf', 'DisplayStyle', 'stairs')
title('hbar_1(151T_C) over 1000 realizations vs Rayleigh pdf')
hold on
a = 0:0.01:3;
plot(a, 2.*a.*exp(-a.^2), 'LineWidth', 1.5);  % Theoretical PDF (page 308, BC)
hold off
legend('hbar_1', 'Rayleigh pdf');
xlabel('a');
ylabel('p_{hbar_1(151T_C)}(a)');