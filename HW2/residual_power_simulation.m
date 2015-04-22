% This script is just a double-check on the theoretical SNR related to the
% truncation of the IR to Nh samples.

%% Clear, initialize useful quantities
% clear all
% close all
% clc

Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Kdb = 3; % 3 dB
K = 10^(Kdb/10);
fd = 5*10^-3/T; % doppler frequency
% these are the coefficients proposed by the authors of anastchugg97
a_dopp = [1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, ...
    -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_dopp = [1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, ...
    6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, ...
    1.8074e-5, 3.0124e-6];
[h_dopp, t_dopp] = impz(b_dopp, a_dopp);
% figure, plot(t_dopp, h_dopp), title('Impulse response of IIR filter h_{ds}')
hds_nrg = sum(h_dopp.^2)
b_dopp = b_dopp / sqrt(hds_nrg);

%% Display the PDP for the channel

% The PDP is the sampling of a continuous time exponential PDP with
% tau_rms / T = 0.3, so tau_rms / Tc = 1.2
% See 4.224 for reference
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time

N_h_max = 900;
residual_power = exp(-N_h_max*Tc/tau_rms);
tau = 0:Tc:N_h_max-1;
pdp_gauss = 1/tau_rms * exp(-tau/tau_rms);

C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|gtilde_i|^2]) = 1 - C^2
pdp_gauss = pdp_gauss.*(1-C^2)/sum(pdp_gauss);

M_d = sum(pdp_gauss);

% Determine suitable length for h, N_h, define criterion
% N_h is final_tau, the criterion is linked to the truncation of
% exp(t/taurms)

pdp = pdp_gauss;
pdp(1) = pdp(1) + C^2;

% Plot the normalised PDP
% figure, stem(tau,10*log10(pdp)), title('Gaussian part of PDP'), xlabel('tau'), ylabel('E[|gtilde_i|^2]');

%% We now have to generate the time behavior of the i-th coefficient. For
% this, we use the model in figure 4.36, page 315.
% The transient will have length equal to the length of the convolution of
% the filter cascade. We will drop that quantity of samples.
% The book suggests using 1/10 < f_d*Tp < 1/5. (Benvenuto said Tp/Tq = 50 or
% more)

Tq = Tc; % Fundamental sampling time. This is the same as Tc, right???
% We stick to anastasopolous chugg paper (1997) and choose Tp such that
% fd*Tp = 0.1
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise
fprintf('fd * Tp = %d \n', Tp*fd);
g_samples_needed = 200000; % Some will be dropped because of transient
w_samples_needed = ceil(g_samples_needed / Tp);
transient = ceil(g_samples_needed/4);

% Doppler filter shape
% [Hf, f] = freqz(b_dopp, a_dopp, 1000, 'whole');
% figure, plot(f/(2*pi), 20*log10(abs(Hf)))

% Generate complex-valued Gaussian white noise with zero mean and unit
% variance
%rng('default');

g_mat = zeros(N_h_max, g_samples_needed - transient);
for ray = 1:N_h_max
    
    w = wgn(w_samples_needed,1,0,'complex');
    % Filter the wgn with a narrowband filter. The filter will have the
    % classical Doppler spectrum in the frequency domain, with f_d * Tc = 1.25*10^-3
    % By using the approach suggested in anastchugg97 we use an iir filter
    % with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
    % Note that it is
    % possible to initialize properly the filters, consider doing it in
    % following versions since it allows to start in steady state conditions
    % and avoid dropping about many samples (it is a very cool thing!)
    gprime = filter(b_dopp, a_dopp, w);
    
    % Interpolation
    t = 1:length(gprime);
    t_fine = Tq/Tp:Tq/Tp:length(gprime);
    
    g_fine = interp1(t, gprime, t_fine, 'spline');
    % Drop the transient
    g_mat(ray, :) = g_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h_max
    g_mat(k, :) = g_mat(k, :)*sqrt(pdp_gauss(k));
end

% Only for LOS component
g_mat(1, :) = g_mat(1, :) + C;

% G = fft(g_mat, [], 2);
% figure, plot(20*log10(abs(G.')));

%% Compute the energy loss if we use a smaller N_h

snr_db = 10;
snr_lin = 10^(snr_db/10);

MAX_NH = 25;
nrg_avg = zeros(MAX_NH, 1);
nrg_across_time = zeros(MAX_NH, size(g_mat, 2));
for N_h = 1:MAX_NH
    % Sum of the powers |g_i(t)|^2 for i=0...Nh-1 (it varies across time)
    nrg_across_time(N_h, :) = sum(abs(g_mat(1:N_h, :)).^2, 1);
    % Average in time of the thing above
    nrg_avg(N_h) = sum(nrg_across_time(N_h, :)) / length(nrg_across_time(N_h, :));
end

% figure, plot(:MAX_NH, nrg_avg), title('Average energy of the IR varying N_h')

% figure, plot(nrg_across_time(2:7, :)')
% title('Energy across time')
% legend('N_h = 2', 'N_h = 3', 'N_h = 4', 'N_h = 5', 'N_h = 6', 'N_h = 7')

% E(|delta h|^2) is the avg energy of the tail we dropped (i.e. of the g_i's for i > Nh-1)
% Plot of residual power
nrg_tail = nrg_avg(MAX_NH) - nrg_avg;
figure, plot(1:MAX_NH, -10*log10(nrg_tail(1:MAX_NH))), grid on

% Plot of Lambda_n
lambda_n = 1./(snr_lin * nrg_tail);
figure, plot(1:MAX_NH, 10*log10(lambda_n)), grid on
title('\Lambda_n')
xlabel('N_h'), ylabel('\Lambda_n [dB]')
ylim([-5 15])


%% PLOT: compare theoretical and simulation results

load lambda_n.mat
figure, hold on
plot(1:length(lambda_n_theoretical), 10*log10(lambda_n_theoretical), 'd-')
plot(1:length(lambda_n_simulation), 10*log10(lambda_n_simulation), 'x--') % We need to cavate this one!!!! Throw it away!!!!!1111
grid on, title('\Lambda_n')
xlabel('N_h'), ylabel('\Lambda_n [dB]')
axis([1 4 -2 10]), ax = gca; ax.XTick = 1:5;
legend('Theoretical', 'Simulation', 'Location', 'NorthWest')