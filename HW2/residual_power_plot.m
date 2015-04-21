clear all;
close all;
clc;

Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
Kdb = 3; % 3 dB
K = 10^(Kdb/10);
C = sqrt(K/(K + 1));
M_d = 1-C^2;
snr = 10; %dB
snr_lin = 10^(snr/10);

% wrong
for N_h = 1:20; % To be determined
    % note that for this choice of final_tau
    h = 1/tau_rms*exp(-(0:899)*Tc/tau_rms);
    h = h/sum(h); % complete & normalized pdp, without truncation
    hhat = zeros(1, 900);
    hhat(1:N_h) = 1/tau_rms*exp(-(0:N_h-1)*Tc/tau_rms);
    hhat = hhat/sum(hhat); % truncated & normalized pdp
    deltah = h - hhat; %NOTE: h = E(|g_i|^2)
    
    power_h(N_h) = sum(abs(h).^2);
    residual_power(N_h) = sum(abs(deltah).^2);
    
    lambda_n(N_h) = 1/(snr_lin * residual_power(N_h));
end

figure, plot(10*log10(power_h./residual_power)), grid on
title('\Lambda_e = |pdp intero|^2/|pdp intero - troncato|^2')
xlabel('N_h')
ylabel('\Lambda_e [dB]')
ylim([0, 20])

figure, plot(10*log10(lambda_n)), grid on
title('\Lambda_n = |1|/(snr*|pdp intero - troncato|^2)')
xlabel('N_h')
ylabel('\Lambda_n [dB]')
ylim([-5, 5])

%%
% What if we consider the sqrt of the exp, since it is the pdp? Note that
% these considerations come from filtering theory, where everything is
% referred to filter coefficients and not to powers

for N_h = 1:20; % To be determined
    % note that for this choice of final_tau
    h = 1/tau_rms*exp(-(0:899)*Tc/tau_rms);
    hhat = zeros(1, 900);
    hhat(1:N_h) = 1/tau_rms*exp(-(0:N_h-1)*Tc/tau_rms);
    deltah = h - hhat; %NOTE: h = E(|g_i|^2)
    
    power_h(N_h) = sum(abs(h));
    residual_power(N_h) = sum(abs(deltah));
    
    lambda_n(N_h) = 1/(snr_lin * residual_power(N_h));
end

figure, plot(10*log10(power_h./residual_power)), grid on
title('\Lambda_e = |sqrt(exp)|^2/|sqrt(exp) - troncato|^2')
xlabel('N_h')
ylabel('\Lambda_e [dB]')
ylim([0, 20])

figure, plot(10*log10(lambda_n)), grid on
title('\Lambda_n = |1|/(snr*|sqrt(exp) - troncato|^2)')
xlabel('N_h')
ylabel('\Lambda_n [dB]')
ylim([-5, 10])

%% Third version
% Each time the exp is truncated it bring a different normalization, thus a
% different error

for N_h = 1:20
    tau = 0:Tc:N_h-1;
    tau_long = 0:Tc:899; % for bigger values MATLAB put the thing to 0
    pdp_gauss = zeros(1, 900);
    pdp_gauss(1:N_h) = 1/tau_rms * exp(-tau/tau_rms);
    pdp_complete = 1/tau_rms * exp(-tau_long/tau_rms);
    
    C = sqrt(K/(K+1));
    % normalize pdp: it must be sum(E[|gtilde_i|^2]) = 1 - C^2
    pdp_gauss = pdp_gauss.*(1-C^2)/sum(pdp_gauss);
    pdp_complete = pdp_complete.*(1-C^2)/sum(pdp_complete);
    
    % compare... what?
    deltapdp = pdp_complete - pdp_gauss;
    residual_power = sum(abs(deltapdp)); % CHECK THIS, no ^2 because already a power??
    lambdapdp_n(N_h) = 1/(snr_lin * residual_power);
end

figure, plot(10*log10(lambdapdp_n)), grid on
title('\Lambda_n = |1|/(snr*|diff pdp|^2)')
xlabel('N_h')
ylabel('\Lambda_n [dB]')
ylim([-5, 5])
