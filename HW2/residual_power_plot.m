clear all;
close all;
clc;

Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time

for N_h = 1:20; % To be determined
    % note that for this choice of final_tau
    residual_power(N_h) = sum(exp(-(N_h:4000)*Tc/tau_rms));
end

figure, plot(10*log10(residual_power))