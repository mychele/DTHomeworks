% This script initializes all system specifications, that will always
% remain constant throughout the execution of all the scripts.

% Sampling times
Tc = 1;     % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Tq = Tc;    % Fundamental sampling time. This is the same as Tc
fd = 5*10^-3/T; % doppler spread
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise, 
% in order to apply Anastasopoulos and Chugg (1997) filter it must be Tp =
% 0.1

Kdb = 3; % 3 dB
K = 10^(Kdb/10);
C = sqrt(K/(K+1)); % this holds if PDP sums to 1

tau_rms = 0.3*T; 

snr = 10; % db
snr_lin = 10^(snr/10);

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