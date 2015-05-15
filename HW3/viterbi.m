%
% q(12:end) = q(12:end) / 20;
% q(1:10) = q(1:10) / 8;

clear
close all
%clc

Tc = 1;
T = 4 * Tc;
L_data = 2^18 - 1;
snr = 6; %dB

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
MEMORY = 20 * N;

%throwitaway:
%newStateBase=@(currstate) (mod(currstate-1, M^(N1+N2-1)) * 4);
%r = [(1 + 1i) * ones(N2, 1); r; (1+1i)*ones(N1,1)]; % add some stuff

% figure
% stem(real(r)), hold on, stem(real(packet), 'x')
% legend('r', 'sent')
% figure, stem(0:N1+N2, real(hi)), hold on, stem(0:N1+N2, imag(hi)), title('hi')

% Init stuff
survSeq = zeros(Ns, min(TMAX, 2*MEMORY));
survSeq(:, 1) = 1:Ns;
survSeq_writingcol = 1;
survSeq_shift = 0;
detectedStates = zeros(1, length(packet));
cost = zeros(Ns, 1); % Define Gamma(-1), i.e. the cost, for each state
%cost = ones(Ns, 1) * Inf;

tic;

for k = 1 : length(r)   % Main loop
    survSeq_writingcol = survSeq_writingcol + 1;
    costnew = - ones(Ns, 1);
    pred = zeros(Ns, 1);
    newstate = 0;
    
    for state = 1 : Ns  % Cycle through all states, at time k(?)
        for j = 1:M     % M possibilities for the new symbol
            
            % Index of the new state: it's mod(state-1, M^(N1+N2-1)) * M + j
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            
            % TODO optimize
            supposednewseq = [symb(mod(survSeq(state, 1:survSeq_writingcol-1)-1,M)+1), symb(j)];
            
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
    % the number of iterations of the inner loops above.
    % TODO optimize the way matrices are handled. Avoid creating so many.
    
    if survSeq_writingcol == size(survSeq, 2)
        % We would write in the last column... shift first!
        
        % Take the first row of old decided states and store it before erasing it.
        detectedStates(1+survSeq_shift : survSeq_shift+MEMORY) = survSeq(1, 1:MEMORY);
        
        % Shift of half the size, that is our memory.
        survSeq(:, 1:MEMORY) = survSeq(:, MEMORY+1:end);
        survSeq_shift = survSeq_shift + MEMORY;
        survSeq_writingcol = MEMORY; % that is 2MEMORY - MEMORY
    end
        
    
    temp = zeros(size(survSeq));
    for newstate = 1:Ns
        temp(newstate, 1:survSeq_writingcol) = ...
            [survSeq(pred(newstate), 1:survSeq_writingcol-1), newstate];
    end
    survSeq = temp;
    
    %allcosts(:, k) = cost;
    cost = costnew;
end

toc
elapsed_time = toc;

% Finish storing, then get the symbols
detectedStates(1+survSeq_shift : survSeq_shift+survSeq_writingcol-1) = ...
    survSeq(1, 1:survSeq_writingcol-1);
detected = symb(mod(detectedStates-1, M) + 1);

num_errors = sum(packet-detected.' ~= 0);
fprintf('P_err is approximately %.g (%d errors)\n', num_errors / length(packet), num_errors)
% TODO: understand why the last symbol is always in error.