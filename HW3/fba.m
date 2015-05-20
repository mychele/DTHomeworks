% TODO functionise me

%function [ detected, pbit ] = fba( packet, r, hi, L1, L2 ) %#ok<INUSL,STOUT>
clear
receiver_util
close all
clc
r = x(1+assumed_dly:1+assumed_dly+length(packet)-1);
M = 4;
L1 = 0;
L2 = 4;
% Fucking ugly, I know.

% This function executes the Max-Log-MAP Algorithm.
% packet: The originally sent data
% r: The received data
% hi: The estimated channel IR
% L1: The number of precursors for each symbol
% L2: The number of postcursors for each symbol

% State definition: s_k = (a_k+L1 ... a_k ... a_k-L2+1)
% States are identified by their indices.

% The desired signal u_k can be obtained as linear combination of s_k
% and s_k-1
% tic
% Initialization
Ns = M^(L1+L2); % Number of states
K = length(r);
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % Possible transmitted symbols (QPSK)

% -- Define u_mat matrix
states_symbols = zeros(Ns, M);
statelength = L1 + L2; % number of digits of the state, i.e. its length in base M=4
statevec = zeros(1, statelength); % symbol index, from the oldest to the newest
u_mat = zeros(Ns, M);
for state = 1:Ns
    % Set value of the current element of u_mat
    for j = 1:M
        lastsymbols = [symb(statevec + 1), symb(j)]; % symbols, from the oldest to the newest
        u_mat(state, j) = lastsymbols * flipud(hi);
    end    
    states_symbols(state,:) = lastsymbols(1:4);
    % Update statevec
    statevec(statelength) = statevec(statelength) + 1;
    i = statelength;
    while (statevec(i) >= M && i > 1)
        statevec(i) = 0;
        i = i-1;
        statevec(i) = statevec(i) + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel transition metrics computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO options: operate on windows, therefore c can be handle and
% precomputed OR try to define it ONLY for the states for which there's a
% transition

c = -1000*ones(Ns, Ns, K+1);
fprintf('channel transition metric...')
for k = 1:K
    newstate = 0;
    
    for state = 1:Ns
        % Iterate over the branches
        for j = 1:M
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            
            % Compute cost
            u_k = u_mat(state, j);
            c(newstate, state, k) = - abs(r(k) - u_k)^2;
        end
    end
end
c(:,:,K+1) = 0;
fprintf('done\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward metric computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(Ns, K+1);   % This will also initialize the last state
fprintf('bck...')
for k = K:-1:1      
    for i = 1:Ns
        b(i, k) = max(b(1:Ns, k+1) + c(1:Ns, i, k+1));
    end
end
fprintf('done\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward metric, state metric, log-like func computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_old = zeros(Ns, 1);   % f_old represents the forward metric at time k-1
f_new = zeros(Ns, 1);   % f_new represents the fwd metric at time k
% v = zeros(Ns, K); % since the decision is taken upon the likelyhood we
% don't need this
l = zeros(M, 1);
decisions = zeros(K, 1);
fprintf('fwd...')
% initialize k = -1!
for j = 1:Ns
    % Only keep the maximum among the fwd metrics
    f_old(j) = max(0 + c(j, 1:Ns, 1)); % 0 until we don't have smth better
    
end
for k = 1:K   % F_(-1) is the initial condition!
    for j = 1:Ns       
        % Only keep the maximum among the fwd metrics
        f_new(j) = max(f_old(1:Ns) + c(j, 1:Ns, k).');
    end
    v = f_new + b(:, k);
    for beta = 1:M
        ind = find(states_symbols(:,M) == symb(beta));
        l(beta) = max(v(ind));
    end
    [~, maxind] = max(l);
    decisions(k) = symb(maxind);
    f_old = f_new;
end
fprintf('done\n')

% disp(toc)

% Pbit computation
[pbit, num_bit_errors] = BER(packet, decisions);

num_errors = sum(packet-decisions ~= 0);
fprintf('P_err = %.g (%d errors)\n', num_errors / length(packet), num_errors)
fprintf('P_bit = %.g (%d errors)\n', pbit, num_bit_errors)

%end