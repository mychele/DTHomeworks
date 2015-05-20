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

c = zeros(M, Ns, K+1);
fprintf('channel transition metric...')
for k = 1:K
    c(:, :, k) = (-abs(r(k) - u_mat).^2).';
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
        % Index of the new state: it's mod(state-1, M^(L1+L2-1)) * M + j
        first_poss_state = mod(i-1, M^(L1 + L2 - 1))*M + 1;
        b(i, k) = max(b(first_poss_state:first_poss_state+M-1, k+1) + c(:, i, k+1));
    end
end
fprintf('done\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward metric, state metric, log-like func computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_old = zeros(Ns, 1);   % f_old represents the forward metric at time k-1, it also initializes k = -1
f_new = zeros(Ns, 1);   % f_new represents the fwd metric at time k
% v = zeros(Ns, K); % since the decision is taken upon the likelyhood we
% don't need this
l = zeros(M, 1);
decisions = zeros(K, 1);
row_step = (0:M-1)*M^(L1+L2-1);
fprintf('fwd...')
% for j = 1:Ns
%     % Only keep the maximum among the fwd metrics
%     f_old(j) = 0; %max(0 + c(j, 1:Ns, 1)); % 0 until we don't have smth better
%     
% end
for k = 1:K   % F_(-1) is the initial condition!
    for j = 1:Ns       
        % Only keep the maximum among the fwd metrics
        % The state l from which I can go to j are ceil(j/M) +
        % r*M^(L1+L2-1), r = 0, 1, .... , M-1
        in_vec = ceil(j/M) + row_step;
        f_new(j) = max(f_old(in_vec) + c(mod(j-1, 4)+1, in_vec, k).');
    end
    v = f_new + b(:, k);
    for beta = 1:M
        ind = find(states_symbols(:,M) == symb(beta));
        l(beta) = max(v(ind));
    end
    [~, maxind] = max(l);
    decisions(k) = symb(maxind); % TODO, check if moving out some things that
    % can be computed with matrices improves the performances
    f_old = f_new;
end
fprintf('done\n')



% toc

% Pbit computation
[pbit, num_bit_errors] = BER(packet, decisions);

num_errors = sum(packet-decisions ~= 0);
fprintf('P_err = %.g (%d errors)\n', num_errors / length(packet), num_errors)
fprintf('P_bit = %.g (%d errors)\n', pbit, num_bit_errors)

%end