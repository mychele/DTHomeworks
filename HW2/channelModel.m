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

% Plot the PDP
figure, stem(tau,pdp), title('PDP'), xlabel('tau'), ylabel('E[|g_i|^2]'), hold on;

C = sqrt(K*sum(pdp)/(K + 1));    % This is the deterministic component inside of g_1
M_d = sum(pdp) - C^2;

%% TODO: Determine suitable length for h, N_h, define criterion

%% Normalize the PDP
% It must be sum(E[|g_i|^2] = 1
pdp = pdp/sum(pdp);

% Plot the normalised PDP
figure, stem(tau,pdp), title('PDP'), xlabel('tau'), ylabel('E[|g_i|^2]');
