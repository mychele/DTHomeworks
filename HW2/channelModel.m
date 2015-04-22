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
[h_dopp, t_dopp] = impz(b_dopp, a_dopp);
figure, plot(t_dopp, h_dopp), title('Impulse response of IIR filter h_{ds}')
hds_nrg = sum(h_dopp.^2)
b_dopp = b_dopp / sqrt(hds_nrg);

%% Display the PDP for the channel

% The PDP is the sampling of a continuous time exponential PDP with
% tau_rms / T = 0.3, so tau_rms / Tc = 1.2
% See 4.224 for reference
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time

N_h = 3; % To be determined
% note that for this choice of final_tau
residual_power = exp(-N_h*Tc/tau_rms);
if (residual_power > 0.1)
    disp('Consider increasing final_tau, since the exp is truncated too early')
end
tau = 0:Tc:N_h-1;
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
figure, stem(tau,10*log10(pdp)), title('Gaussian part of PDP'), xlabel('tau'), ylabel('E[|gtilde_i|^2]');

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
[Hf, f] = freqz(b_dopp, a_dopp, 1000, 'whole');
figure, plot(f/(2*pi), 20*log10(abs(Hf)))
title('Frequency response of the Doppler filter')

% Generate complex-valued Gaussian white noise with zero mean and unit
% variance
rng('default');

g_mat = zeros(N_h, g_samples_needed - transient);
for ray = 1:N_h
    
    w = wgn(w_samples_needed,1,0,'complex');
    fprintf('variance of white noise=%d \n', var(w))
    % Filter the wgn with a narrowband filter. The filter will have the
    % classical Doppler spectrum in the frequency domain, with f_d * Tc = 1.25*10^-3
    % By using the approach suggested in anastchugg97 we use an iir filter
    % with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
    % Note that it is
    % possible to initialize properly the filters, consider doing it in
    % following versions since it allows to start in steady state conditions
    % and avoid dropping about many samples (it is a very cool thing!)
    gprime = filter(b_dopp, a_dopp, w);
    fprintf('g %d after iterpolation has mean %d and variance %d  \n', ray, mean(gprime), var(gprime));
    %Gpr = fft(gprime);
    %figure, plot(20*log10(abs(Gpr)))
    
    % Interpolation
    t = 1:length(gprime);
    t_fine = Tq/Tp:Tq/Tp:length(gprime);
    
    g_fine = interp1(t, gprime, t_fine, 'spline');
    fprintf('g %d after iterpolation has mean %d and variance %d  \n', ray, mean(g_fine), var(g_fine));
    % Drop the transient
    g_mat(ray, :) = g_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h
    g_mat(k, :) = g_mat(k, :)*sqrt(pdp_gauss(k));
end

% Only for LOS component
g_mat(1, :) = g_mat(1, :) + C;

G = fft(g_mat, [], 2);

figure, plot((0:size(G, 2)-1)/size(G,2), 20*log10(abs(G.')));
title('Magnitude of DFT of g_i')

%% Show the behavior of |g_1(nTc)| for n = 0:1999, dropping the transient

figure, hold on
plot(0:1999, 20*log10(abs(g_mat(2, 1:2000).')))
grid on, box on, xlabel('Time samples'), ylabel('|g_1(nT_C)|_{dB}')
title('|g_1(nT_C)|')


%% Show all |g_i|'s

figure, hold on
plot(0:9999, 20*log10(abs(g_mat(:, 1:10000).')))
grid on, box on, xlabel('Time samples'), ylabel('|g_i(nT_C)|_{dB}')
legend('g_0', 'g_1', 'g_2'), ylim([-40 10])
title('|g_i|')

% g_mean = mean(g_mat(:, 1:5000).');
% g_var = 10*log10(var(g_mat(:, :).'));



%% Histogram of g_1

figure, histogram(abs(g_mat(2, 1:1000).')/sqrt(pdp_gauss(2)), 100)
title('1000 samples of |g_1|')

%% Histograms of g_1 in subsequent disjoint intervals
% The variability shows that these intervals of 1000 samples are too small.
% Looking at |g1| in one of these intervals, we observe that this quantity
% has an oscillatory trend. In fact, g1(nTc) for different n's are
% identically distributed rv's, but they are not independent. Whenever
% there is a peak (max or min), the values of g1 for consecutive
% time instants are similar, therefore they build up on the histogram,
% showing a small number of peaks in the histogram itself. These peaks are
% the values of g1 around the instants nTc in which it has a local max or min.
% If we take a sufficiently large time window, the number of maxima and
% minima increases, and they can be considered independent if they are
% enough (?).
% Considering two different realizations of the channel impulse response,
% two samples of g1 at the same time (taken from these two realizations)
% are not only identically distributed, but also independent. Therefore
% with as little as 1000 such samples we observe a distribution very close
% to a Rayleigh distribution, as we expect from our model.

% Uncomment the following code to see that variability

% figure
% width = 1000;
% for i = 0:9
%     histogram(abs(g_mat(2, (i*width + 1):((i+1)*width)).')/sqrt(pdp_gauss(2)), 100)
%     axis([0 2.5 0 0.05*width])
%     pause
% end

%% Autocorrelation of g_1

begin = 1;
end_factor = 0.2;
g1_temp = g_mat(2, begin:end*end_factor).';
N_corr = 4000; %floor(length(g1_temp/10)); % there are lot of samples!

autoc_1_intime = autocorrelation(g1_temp, N_corr);
z = (0:length(autoc_1_intime)-1)*2*pi*fd; %support for bessel function, LAG not actual time
M_1 = autoc_1_intime(1); %stat power of g_1 over the considered support in time
bess = besselj(0, z); %bessel function of order 0

figure, plot(real(autoc_1_intime)), hold on, plot(bess*M_1)
legend('Real(r_g(\Delta t))', 'Bessel function, first type, order 0')
title('Autoc from t=0')

begin = 30000;
end_factor = 0.5;
g1_temp = g_mat(2, begin:end*end_factor).';
N_corr = 4000; %floor(length(g1_temp/10)); % there are lot of samples!

autoc_1_intime = autocorrelation(g1_temp, N_corr);
z = (0:length(autoc_1_intime)-1)*2*pi*fd;
M_1 = autoc_1_intime(1);
bess = besselj(0, z);

figure, plot(real(autoc_1_intime)), hold on, plot(bess*M_1)
legend('Real(r_g(\Delta t))', 'Bessel function, first type, order 0')
title(sprintf('Autoc from t=%d', begin))


%% Simulation

rng('default') % for reproducibility
numexp = 1000;
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
title('|g1(151T_C)| over 1000 realizations')