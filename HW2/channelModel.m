%% Clear, initialize useful quantities
clear all
close all
clc

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

%% Display the PDP for the channel

% The PDP is the sampling of a continuous time exponential PDP with
% tau_rms / T = 0.3, so tau_rms / Tc = 1.2
% See 4.224 for reference
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time

N_h = 5; % To be determined
% note that for this choice of final_tau
residual_power = exp(-N_h*Tc/tau_rms);
if (residual_power > 0.1)
    disp('Consider increasing final_tau, since the exp is truncated too early')
end
tau = 0:Tc:N_h;
pdp = 1/tau_rms * exp(-tau/tau_rms);

% normalize pdp: it must be sum(E[|g_i|^2] = 1
pdp = pdp/sum(pdp);

C = sqrt(K/(K + 1));   % This is the deterministic component inside of g_1
M_d = sum(pdp) - C^2; % 1 - C^2, it's the sum of energies of random components gtilde

% Determine suitable length for h, N_h, define criterion
% N_h is final_tau, the criterion is linked to the truncation of
% exp(t/taurms)

% Plot the normalised PDP
figure, stem(tau,10*log10(pdp)), title('PDP'), xlabel('tau'), ylabel('E[|g_i|^2]');

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

% Generate complex-valued Gaussian white noise with zero mean and unit
% variance
rng('default');
w = wgn(w_samples_needed,1,0,'complex');

% Filter the wgn with a narrowband filter. The filter will have the
% classical Doppler spectrum in the frequency domain, with f_d * Tc = 1.25*10^-3
% By using the approach suggested in anastchugg97 we use an iir filter
% with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
% Note that it is
% possible to initialize properly the filters, consider doing it in
% following versions since it allows to start in steady state conditions
% and avoid dropping about 2000 samples (it is a very cool thing!)

[Hf, f] = freqz(b_dopp, a_dopp, 1000, 'whole');
figure, plot(f/(2*pi), 20*log10(abs(Hf)))

gprime = filter(b_dopp, a_dopp, w);

Gpr = fft(gprime);
figure, plot(20*log10(abs(Gpr)))

% Interpolation
t = 1:length(gprime);
t_fine = Tq/Tp:Tq/Tp:length(gprime);

g_fine = interp1(t, gprime, t_fine, 'spline');

% Drop the transient
g_notrans = g_fine(50000:end);

% Energy scaling
g_mat = repmat(g_notrans, N_h + 1, 1);

for k = 1:N_h
    g_mat(k, :) = g_mat(k, :)*sqrt(pdp(k));
end

% Only for LOS component
g_mat(1, :) = g_mat(1, :) + C;

G = fft(g_mat, [], 2);

figure, plot(20*log10(abs(G.')));

%% Show the behavior of |g_1(nTc)| for n = 0:1999, dropping the transient

figure, plot(0:1999, abs(g_mat(2, 1:2000).'))
g_mean = mean(g_mat(:, 1:5000).');

%% Histogram of g_1

figure, histogram(abs(g_mat(2, 1:1000).')/sqrt(pdp(2)), 100)

%% Simulation

rng('default') % for reproducibility
numexp = 5000;
g_1 = zeros(numexp, 1);
for k = 1:numexp
    disp(k)
    w = wgn(w_samples_needed,1,0,'complex');
        
    gprime = filter(b_dopp, a_dopp, w);
    
    % Interpolation
    t = 1:length(gprime);
    t_fine = Tq/Tp:Tq/Tp:length(gprime);
    
    g_fine = interp1(t, gprime, t_fine, 'spline');
    
    % Drop the transient
    g_notrans = g_fine(50000:end);
    
    % Energy scaling
    g_1(k) = g_notrans(152);  % it should be multiplied by sqrt(pdp(2)) but since
    % for the histogram it is required to divide for the same factor then
    % we drop this computation to save computational time
end

figure, histogram(abs(g_1), 50)


