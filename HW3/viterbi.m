%
% q(12:end) = q(12:end) / 20;
% q(1:10) = q(1:10) / 8;

clear
close all
%clc

Tc = 1;
T = 4 * Tc;
L_data = 2^20 - 1;
snr = 50; %dB

% create, send and receive data with the given channel
fprintf('\nCalling txrc()\n')
fprintf('Generating input symbols and channel output... ')
[packet, r_T4, ~] = txrc(L_data, snr, T, Tc);
fprintf('done!\n')

% estimate the channel using the first 100 samples (4*length(ts))
N = 3;
fprintf('Calling get_channel_info()\n')
fprintf('Estimating timing phase and IR... ')
[ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);
fprintf('done!\n')

% sample to get r @ T
init_offs = mod(m_opt, 4); % in T/4
t0 = floor(m_opt / 4); % from now consider T = 1
r = r_T4(init_offs+1:T:end); % data sampled in T
r = r / h(N1+1); % data normalized by h0
r = r(N1+1 : end-N2);    %discard precursors and postcursors (initial and final samples)
hi = h/h(N1+1); % impulse response normalized by h0

%% Ready. Do Viterbi.

% State: most recent is at the left, as in the book, and it has lowest
% weight.

fprintf('\nViterbi started.\n\n')

M = 4;
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % possible transmitted symbols (QPSK)
Ns = M ^ (N1+N2);
TMAX = L_data + 50;

%throwitaway:
%newStateBase=@(currstate) (mod(currstate-1, M^(N1+N2-1)) * 4);
%r = [(1 + 1i) * ones(N2, 1); r; (1+1i)*ones(N1,1)]; % add some stuff

% figure
% stem(real(r)), hold on, stem(real(packet), 'x')
% legend('r', 'sent')
% figure, stem(0:N1+N2, real(hi)), hold on, stem(0:N1+N2, imag(hi)), title('hi')

% Init stuff
surv_seq = zeros(Ns, TMAX);
surv_seq(:, 1) = 1:Ns;
surv_seq_len = 1;
cost = zeros(Ns, 1); % Define Gamma(-1), i.e. the cost, for each state
%cost = ones(Ns, 1) * Inf;

tic;

for k = 1 : length(r)   % Main loop
    surv_seq_len = surv_seq_len + 1;
    costnew = - ones(Ns, 1);
    pred = zeros(Ns, 1);
    newstate = 0;
    k
    
    for state = 1 : Ns  % Cycle through all states, at time k(?)
        for j = 1:M     % M possibilities for the new symbol
            
            % Index of the new state: it's mod(state-1, M^(N1+N2-1)) * M + j
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            
            % TODO optimize
            supposednewseq = [symb(mod(surv_seq(state, 1:surv_seq_len-1)-1,M)+1), symb(j)];
            
            % TODO optimize
            difflength = N - length(supposednewseq);
            if difflength > 0
                supposednewseq = [zeros(1, difflength), supposednewseq];
            else
                supposednewseq = supposednewseq(end-N+1:end);
            end
            u = supposednewseq * flipud(hi);
            newstate_cost = cost(state) + abs(r(k) - u)^2;
            if costnew(newstate) == -1 ...     % not assigned yet, or...
                    || costnew(newstate) > newstate_cost    % ...found path with lower cost
                % Update everything
                costnew(newstate) = newstate_cost;
                pred(newstate) = state;
            end
            
%             fprintf('k=%d\t%d->%d\tu=%f+j%f\tr=%f+j%f\n', k, state, newstate, real(u), imag(u), real(r(k)), imag(r(k)))
            
        end
    end
    
    
    % The following operations strongly affect the computation time, if
    % the number of states is not too large. Otherwise the bottleneck is
    % the number of iterations of the inner loops above. This stuff is slow
    % with long sequences since it handles huge matrices. Split them!
    % TODO optimize the way matrices are handled. Avoid creating so many!
    
    temp = zeros(size(surv_seq));
    for newstate = 1:Ns
        temp(newstate, 1:surv_seq_len) = [surv_seq(pred(newstate), 1:surv_seq_len-1), newstate];
    end
    surv_seq = temp;
    
    allcosts(:, k) = cost;
    cost = costnew;
    temp = zeros(size(surv_seq));
    for l = 1 : Ns
        index = surv_seq(l, surv_seq_len);
        if index > 0
            temp(surv_seq(l, surv_seq_len), :) = surv_seq(l, :);
        end
    end
    surv_seq = temp;
    
end

toc

% Take the first row and get the symbols
detected = symb(mod(surv_seq(1, 1:surv_seq_len-1)-1, M) + 1);

fprintf('Pbit is approximately %.g\n', sum(packet-detected.' ~= 0) / length(packet))