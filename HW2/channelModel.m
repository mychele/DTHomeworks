%% Clear, initialize useful quantities
clear all
close all
clc

Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
K = 2; % 3 dB

%% Display the PDP for the channel

% The PDP is the sampling of a continuous time exponential PDP with 
% tau_rms / T = 0.3, so tau_rms / Tc = 1.2
% See 4.224 for reference
final_tau = 5; % To be determined
tau = 0:Tc:final_tau;
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
pdp = 1/tau_rms * exp(-tau/tau_rms);

C = sqrt(K*sum(pdp)/(K + 1));    % This is the deterministic component inside of g_1
M_d = sum(pdp) - C^2;

%% TODO: Determine suitable length for h, N_h, define criterion

%% Normalize the PDP
% It must be sum(E[|g_i|^2] = 1
pdp = pdp/sum(pdp);

% Plot the normalised PDP
figure, stem(tau,10*log10(pdp)), title('PDP'), xlabel('tau'), ylabel('E[|g_i|^2]');

%% Show the behavior of |g_i(nTc)| for n = 0:1999, dropping the transient
% We now have to generate the time behavior of the i-th coefficient. For
% this, we use the model in figure 4.36, page 315.
% The transient will have length equal to the length of the convolution of
% the filter cascade. We will drop that quantity of samples.
% The book suggests using 1/10 < f_d*Tp < 1/5. For f_d*Tp = 1/7.5 we should
% pick Tp = 26.66 T ~= 25 T
Tq = Tc; % Fundamental sampling time. This is the same as Tc, right???
Tp = 25 * Tq; % Sampling time used for filtering the white noise
g_samples_needed = 2000;
w_samples_needed = g_samples_needed / Tp;

% Generate complex-valued Gaussian white noise with zero mean and unit
% variance
w = wgn(w_samples_needed,1,0,'complex');

% Filter the wgn with a narrowband filter. The filter will have the
% classical Doppler spectrum in the frequency domain, with f_d * T = 5*10^-3
