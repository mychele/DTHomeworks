%% Impulse response estimation
% We are at the receiver, we know what the sender is sending and we try to
% estimate it with the LS method (for reference, see page 244).

% Note: As the receiver, we do _not_ know neither N_h nor sigma_w

clear all
close all
clc

% To observe L samples, we need to send L+N-1 samples of the training
% sequence {x(0), ..., x((N-1)+(L-1))}
% TODO Try multiple combinations. Generally, L = 2*N
L = 3; % Length of the observation 
N = 1; % Length of the impulse response of the channel

%% Generate training sequence
% The x sequence must be a partially repeated M-L sequence of length L. We
% need it to have size L+N-1.
% this holds for L = 15
r = log2(L+1);
p = zeros(L,1);
p(1:r) = ones(1,r).'; % Set arbitrary initial condition
for l = r+1:(L)
    p(l) = xor(p(l-1), p(l-2));
end
clear l

x = [p; p(1:N-1)];
x(x == 0) = -1;

% Construct a PN object
% h = commsrc.pn('GenPoly', [4 3 0]);
% set(h, 'NumBitsOut', 1);
% set(h, 'InitialStates', ones(4,1));
% 
% mls = zeros(1,L);
% for k = 1:L
%     mls(k) = generate(h);
% end

%% Generate gi

N_h = 4; 
a_dopp = [1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, ...
    -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_dopp = [1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, ...
    6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, ...
    1.8074e-5, 3.0124e-6];
[h_dopp, t_dopp] = impz(b_dopp, a_dopp);
Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Kdb = 3; % 3 dB
K = 10^(Kdb/10);
fd = 5*10^-3/T; % doppler frequency
Tq = Tc; % Fundamental sampling time. This is the same as Tc, right???
% We stick to anastasopolous chugg paper (1997) and choose Tp such that
% fd*Tp = 0.1
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
tau = 0:Tc:N_h-1;
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise
pdp_gauss = 1/tau_rms * exp(-tau/tau_rms);

C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|gtilde_i|^2]) = 1 - C^2
pdp_gauss = pdp_gauss.*(1-C^2)/sum(pdp_gauss);

M_d = sum(pdp_gauss);

snr = 10; % db
snr_lin = 10^(snr/10);
sigma_w = 1/(T/Tc*snr);


fprintf('fd * Tp = %d \n', Tp*fd);
g_samples_needed = 200000; % Some will be dropped because of transient
w_samples_needed = ceil(g_samples_needed / Tp);
transient = ceil(g_samples_needed/4);

% Generate complex-valued Gaussian white noise with zero mean and unit
% variance
%rng('default');
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

%% Generate desired signal via a polyphase implementation
d = zeros(T*length(x),1);
for k = 0:length(x)-1
    % Generate white noise
    w = wgn(4, 1, 10*log10(1/40), 'complex');
    for idx = 0:N_h-1
        d(k*T + idx*Tc + 1) = g_mat(idx+1,k*T + idx*Tc+1) * x(k+1) + w(idx+1); 
    end
end

figure, stem(0:T:(length(d)-1), abs(x)), hold on, stem(0:Tc:length(d)-1, abs(d));
legend('x', 'd');

% Using the data matrix (page 246), easier implementation
h_hat = zeros(N_h-1,1);
for idx = 0:(N_h-1)
    I = zeros(L,N);
    for column = 1:N
        I(:,column) = x(N-column+1:(N+L-column));
    end
    o = d((idx + 1):(N_h):end);
    
    Phi = I'*I;
    theta = I'*o;
    
    h_hat(idx+1) = inv(Phi) * theta.';
end
figure, stem(abs(g_mat(:,1))), hold on, stem(abs(h_hat));
legend('h', 'h_hat');

% NOTE: the training sequence is defined in T, while the channel and the
% receiver operate in Tc = T/4
% Thus it is necessary to define the polyphase realization of the system
% and perform 4 estimation problems "in parallel"

% A trivial (very bad for Benvenuto) implementation in theory would be
% d0(kT + 0Tc) = g0(kT + 0Tc)x(kT) + eventually other terms at distance 4
% d1(kT + 1Tc) = g1(kT + 1Tc)x(kT)
% d2(kT + 2Tc) = g1(kT + 2Tc)x(kT)
% d3(kT + 1Tc) = g1(kT + 3Tc)x(kT) = 0 always?

% h = rand(N,1);
% d = zeros(N+L, 1);
% w = rand(N+L, 1);
% for k = N:(N+L)
%     %d(k) 
% end


% NITER = 5000;
% for iii=1:NITER
%     for k = 0:length(x)-1
%         % Generate white noise
%         w = wgn(4, 1, 10*log10(1/40), 'complex');
%         
%         % Generate d
%         d(k*T + 1) = g_mat(1,(iii-1)*length(x)*T + k*T+1)*x(k+1) + w(1);
%         d(k*T + Tc + 1) = g_mat(2,(iii-1)*length(x)*T +k*T + Tc+1) * x(k+1) + w(2);
%         d(k*T + 2*Tc + 1) = g_mat(3,(iii-1)*length(x)*T +k*T + 2*Tc+1) * x(k+1) + w(3);
%         d(k*T + 3*Tc + 1) = w(4); %g_mat(4,k*T + 3*Tc+1) * x(k+1) + w(4);
%         
%         % Using the data matrix (page 246), easier implementation
%         I = zeros(L,N);
%         for column = 1:N
%             I(:,column) = x(N-column+1:(N+L-column));
%         end
%         o = d(N:(N+L-1));
%         
%         Phi = I'*I;
%         theta = I'*o;
%         
%         h_hat = inv(Phi) * theta;
%     end
%     davg = davg + abs(d).^2/NITER;
% end
% figure, stem(abs(davg))

% Using the data matrix (page 246), easier implementation
%         I = zeros(L,N);
%         for column = 1:N
%             I(:,column) = x(N-column+1:(N+L-column));
%         end
%         o = d(N:(N+L-1));
%         
%         Phi = I'*I;
%         theta = I'*o;
%         
%         h_hat(idx+1, k*T + idx*Tc +1) = inv(Phi) * theta;        