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
L = 15; % Length of the observation. 
N = 7; % Length of the impulse response of the channel

%% Generate training sequence
% The x sequence must be a partially repeated M-L sequence of length L. We
% need it to have size L+N-1.
% this holds for L = 15
r = log2(L+1);
p = zeros(L,1);
p(1:r) = [0 0 0 1].'; % Set arbitrary initial condition
for l = r+1:(L)
    p(l) = xor(p(l-3), p(l-4));
end
clear l

x = [p; p(1:N-1)];

%% Generate desired signal
% d = h*x + w

% NOTE: the training sequence is defined in T, while the channel and the
% receiver operate in Tc = T/4
% Thus it is necessary to define the polyphase realization of the system
% and perform 4 estimation problems "in parallel"

% A trivial (very bad for Benvenuto) implementation in theory would be
% d0(kT + 0Tc) = g0(kT + 0Tc)x(kT) + eventually other terms at distance 4
% d1(kT + 1Tc) = g1(kT + 1Tc)x(kT)
% d2(kT + 2Tc) = g1(kT + 2Tc)x(kT)
% d3(kT + 1Tc) = g1(kT + 3Tc)x(kT) = 0 always?

h = rand(N,1);
d = zeros(N+L, 1);
w = rand(N+L, 1);
for k = N:(N+L)
    %d(k) 
end

%% Set up the receiver in order to estimate the channel coefficients
% Using the data matrix (page 246), easier implementation
I = zeros(L,N);
for column = 1:N
    I(:,column) = x(N-column+1:(N+L-column));
end
o = d(N:(N+L));
Phi = I'*I;

