%% Data

% Sampling times
Tc = 1;     % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Tq = Tc;    % Fundamental sampling time. This is the same as Tc
fd = 5*10^-3/T; % doppler spread
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise, 
% in order to apply Anastasopoulos and Chugg (1997) filter it must be Tp =
% 0.1

% PDP (aleatory part)
N_h = 3;
Kdb = 3;    % 3 dB
K = 10^(Kdb/10);
tau_rms = 0.3*T; 
tau = 0:Tc:N_h-1;
M_iTc = 1/tau_rms * exp(-tau/tau_rms);
C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|htilde_i|^2]) = 1 - C^2
M_iTc = M_iTc.*(1-C^2)/sum(M_iTc);

% snr and noise power
snr = 10; % db
snr_lin = 10^(snr/10);
sigma_w = 1/(T/Tc*snr); % the PN sequence has power 1

% Filter for doppler spectrum. The filter will have the
% classical Doppler spectrum in the frequency domain, with f_d * T = 5*10^-3
% By using the approach suggested in Anastasopoulos and Chugg (1997) we use an iir filter
% with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
a_dopp = [1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, ...
    -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_dopp = [1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, ...
    6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, ...
    1.8074e-5, 3.0124e-6];
% The energy needs to be normalized to 1
[h_dopp, ~] = impz(b_dopp, a_dopp);
hds_nrg = sum(h_dopp.^2);
b_dopp = b_dopp / sqrt(hds_nrg);

%% Generation of impulse responses
h_samples_needed = 2000000; % Some will be dropped because of transient, since
% enough time, memory and computational power are available 
w_samples_needed = ceil(h_samples_needed / Tp);
% The filter is IIR, from Anastasopoulos and Chugg (1997) it appears that 
% the effect of the transient is present in about 2000 samples for an interpolation
% factor Q = 100. This model uses Q = 80, since memory and computational power
% are not an issue, in order to be conservative it drops 500000 samples. 
transient = ceil(h_samples_needed/4);

h_mat = zeros(N_h, h_samples_needed - transient);
for ray = 1:N_h
    % Generate complex-valued Gaussian white noise with zero mean and unit variance
    w = wgn(w_samples_needed,1,0,'complex');
    %fprintf('variance of white noise=%d \n', var(w))
    
    hprime = filter(b_dopp, a_dopp, w);
    %fprintf('h%d before interpolation has mean %d and variance %d  \n', ray, mean(hprime), var(hprime));
    
    % Interpolation
    t = 1:length(hprime);
    t_fine = Tq/Tp:Tq/Tp:length(hprime);
    
    h_fine = interp1(t, hprime, t_fine, 'spline');
    %fprintf('h%d after interpolation has mean %d and variance %d  \n', ray, mean(h_fine), var(h_fine));
    % Drop the transient
    h_mat(ray, :) = h_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h
    h_mat(k, :) = h_mat(k, :)*sqrt(M_iTc(k));
end

% Only for LOS component
h_mat(1, :) = h_mat(1, :) + C;

clear tau tau_rms a_dopp b_dopp g_fine gprime h_dopp pdp_gauss t_dopp t_fine g_samples_needed w_samples_needed M_d k ray
