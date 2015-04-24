N_h = 3;
a_dopp = [1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, ...
    -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_dopp = [1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, ...
    6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, ...
    1.8074e-5, 3.0124e-6];
[h_dopp, ~] = impz(b_dopp, a_dopp);
hds_nrg = sum(h_dopp.^2);
b_dopp = b_dopp / sqrt(hds_nrg);
Tc = 1;     % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Kdb = 3;    % 3 dB
K = 10^(Kdb/10);
fd = 5*10^-3/T; % doppler frequency
Tq = Tc;    % Fundamental sampling time. This is the same as Tc, right???
% We stick to anastasopolous chugg paper (1997) and choose Tp such that
% fd*Tp = 0.1
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
tau = 0:Tc:N_h-1;
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise
pdp_gauss = 1/tau_rms * exp(-tau/tau_rms);

C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|gtilde_i|^2]) = 1 - C^2
pdp_gauss = pdp_gauss.*(1-C^2)/sum(pdp_gauss);

snr = 10; % db
snr_lin = 10^(snr/10);
sigma_w = 1/(T/Tc*snr); % the PN sequence has power 1


fprintf('fd * Tp = %d \n', Tp*fd);
g_samples_needed = 200000; % Some will be dropped because of transient
w_samples_needed = ceil(g_samples_needed / Tp);
transient = ceil(g_samples_needed/4);

% Generate complex-valued Gaussian white noise with zero mean and unit variance
%rng('default');
g_mat = zeros(N_h, g_samples_needed - transient);
for ray = 1:N_h
    w = wgn(w_samples_needed,1,0,'complex');
    %fprintf('variance of white noise=%d \n', var(w))
    
    % Filter the wgn with a narrowband filter. The filter will have the
    % classical Doppler spectrum in the frequency domain, with f_d * Tc = 1.25*10^-3
    % By using the approach suggested in anastchugg97 we use an iir filter
    % with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
    % Note that it is
    % possible to initialize properly the filters, consider doing it in
    % following versions since it allows to start in steady state conditions
    % and avoid dropping about many samples (it is a very cool thing!)
    gprime = filter(b_dopp, a_dopp, w);
    %fprintf('g%d before interpolation has mean %d and variance %d  \n', ray, mean(gprime), var(gprime));
    
    % Interpolation
    t = 1:length(gprime);
    t_fine = Tq/Tp:Tq/Tp:length(gprime);
    
    g_fine = interp1(t, gprime, t_fine, 'spline');
    %fprintf('g%d after interpolation has mean %d and variance %d  \n', ray, mean(g_fine), var(g_fine));
    % Drop the transient
    g_mat(ray, :) = g_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h
    g_mat(k, :) = g_mat(k, :)*sqrt(pdp_gauss(k));
end

% Only for LOS component
g_mat(1, :) = g_mat(1, :) + C;

clear tau tau_rms a_dopp b_dopp g_fine gprime h_dopp pdp_gauss t_dopp t_fine g_samples_needed w_samples_needed M_d k ray
