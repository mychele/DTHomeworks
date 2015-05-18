%function [ detected, pbit ] = fba( packet, r, hi, L1, L2 ) %#ok<INUSL,STOUT>
clear all
receiver_util 
close all
clc
r = x(1+assumed_dly:end);
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

% Initialization
Ns = M^(L1+L2); % Number of states
K = length(r);
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % Possible transmitted symbols (QPSK)

% -- Define u_mat matrix
ndigits = L1 + L2; % n digits of the state
statevec = zeros(1, ndigits); % symbol index, from the oldest to the newest
i = ndigits;
states_symbols = zeros(Ns, M);
u_mat = zeros(Ns, M);
for state = 1:Ns
    for j = 1:M
        lastsymbols = [symb(statevec + 1), symb(j)];   % symbols, from the oldest to the newest
        u_mat(state, j) = lastsymbols * flipud(hi);
    end
    
    states_symbols(state,:) = lastsymbols(1:4);
    
    % Update statevec
    statevec(i) = statevec(i) + 1;
    l = i;
    while (statevec(l) >= M && l > 1)
        statevec(l) = 0;
        l = l-1;
        statevec(l) = statevec(l) + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel transition metrics computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = -1000*ones(Ns, Ns, K+1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward metric computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(Ns, K+1);   % This will also initialize the last state

for k = K:-1:1
    
    newstate = 0;
    
    for i = 1:Ns
        
        curr_bkw_metric = zeros(1,Ns);
        
        % Iterate over the trellis transitions for the current state
        for m = 1:Ns
            
            % Compute backward metric for each of the new states
            curr_bkw_metric(m) = b(m, k+1) + c(m, i, k+1);
        end
        
        % Only keep the maximum among the bkw metrics
        b(i,k) = max(curr_bkw_metric);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward metric computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(Ns, K+2);   % This will also initialize the initial state condition
for k = 2:K+1   % F_(-1) is the initial condition!
    
    for j = 1:Ns
        
        curr_fwd_metric = zeros(1,Ns);
        
        % Iterate over the branches
        for l = 1:Ns
            % Compute backward metric for each of the new states
            curr_fwd_metric(l) = f(l, k-1) + c(j, l, k);
        end
        
        % Only keep the maximum among the bkw metrics
        f(j,k) = max(curr_fwd_metric);
    end
end

% State metric computation
v = zeros(Ns, K);
for k = 1:K
    for i = 1:Ns
        v(i, k) = f(i, k+1) + b(i, k);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Likelihood function computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = zeros(M, K);
for k = 1:K
    for beta = 1:M
        l(beta, k) = max(v(find(states_symbols(:,1) == symb(beta)),k));
    end
end

% Decision
[maximum, maxind] = max(l);
decisions = symb(maxind);


% Pbit computation
detected = decisions(3:end-2);
[pbit, num_bit_errors] = BER(packet, detected);

num_errors = sum(packet-detected.' ~= 0);
fprintf('P_err = %.g (%d errors)\n', num_errors / length(packet), num_errors)
fprintf('P_bit = %.g (%d errors)\n', pbit, num_bit_errors)

%end