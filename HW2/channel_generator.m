%% Data and variables initialization

data_init;

% PDP (aleatory part)
N_h = 3;
tau = 0:Tc:N_h-1;
M_iTc = 1/tau_rms * exp(-tau/tau_rms);
C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|htilde_i|^2]) = 1 - C^2
M_iTc = M_iTc.*(1-C^2)/sum(M_iTc);

%% Impulse responses generation
% Some will be dropped because of transient, since enough time, memory and 
% computational power are available. 
h_samples_needed = 8000000 + ceil(Tp/Tc*length(h_dopp)); 
w_samples_needed = ceil(h_samples_needed / Tp);
% The filter is IIR, from Anastasopoulos and Chugg (1997) it appears that 
% the effect of the transient is present in about 2000 samples for an 
% interpolation of factor Q = 100. This model uses Q = 80. Since memory and 
% computational power are not an issue, in order to be conservative it 
% drops 80*length(h_dopp) samples. 
transient = ceil(Tp/Tc*length(h_dopp));

h_mat = zeros(N_h, h_samples_needed - transient);
for ray = 1:N_h
    % Complex-valued Gaussian white noise with zero mean and unit variance
    w = wgn(w_samples_needed,1,0,'complex');
    
    hprime = filter(b_dopp, a_dopp, w);
    
    % Interpolation
    t = 1:length(hprime);
    t_fine = Tq/Tp:Tq/Tp:length(hprime);
    h_fine = interp1(t, hprime, t_fine, 'spline');
    
    % Drop the transient
    h_mat(ray, :) = h_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h
    h_mat(k, :) = h_mat(k, :)*sqrt(M_iTc(k));
end

% Only for LOS component, add deterministic component
h_mat(1, :) = h_mat(1, :) + C;

clear tau tau_rms h_fine hprime pdp_gauss t_dopp t_fine ...
    h_samples_needed w_samples_needed M_d k ray
