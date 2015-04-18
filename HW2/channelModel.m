%% Clear, initialize useful quantities
Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
K = 2; % 3 dB

%% Display the PDP for the channel

% The PDP is the sampling of a continuous time exponential PDP with 
% tau_rms / T = 0.3, so tau_rms / Tc = 0.9
% See 4.224 for reference
final_tau = 5; % To be determined
tau = 0:Tc:final_tau;
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
pdp = 1/tau_rms * exp(-tau/tau_rms);

M_d = sum(pdp);     % This is the power of the PDP
C = sqrt(K*M_d);    % This is the deterministic component inside of g_1

%% TODO: Determine suitable length for h, N_h, define criterion

%% Normalize the PDP
% It must be sum(E[|g_i|^2] = 1
pdp = pdp/M_d;

% Plot the PDP
figure, stem(tau,pdp), title('PDP'), xlabel('tau'), ylabel('E[|g_i|^2]');
