function [ detected, pbit, num_bit_errors, Jmin ] = viterbi( packet, r, hi, N1, N2 , t0 )
%VITERBI

% TODO: use t0 to ignore precursors by maintaining delay of the channel.
% Right now the number of precursors considered is actually N1 that is t0.

% TODO: understand why the last symbol is always in error.


% State: most recent is at the left, as in the book, and it has lowest weight
% (thus, in the base-M representation, it is of course the rightmost one).

fprintf('Viterbi started.\n')


% --- Setup

M = 4;
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % possible transmitted symbols (QPSK)
Ns = M ^ (N1+N2);       % Number of states
TMAX = length(r) + 100; % Max value of time k
N = N1 + N2 + 1;
MEMORY = 20 * N;        %



% --- Init stuff

tic;
survSeq = zeros(Ns, min(TMAX, 2*MEMORY));
survSeq(:, 1) = 1:Ns;
survSeq_writingcol = 1;
survSeq_shift = 0;
detectedStates = zeros(1, length(packet));
cost = zeros(Ns, 1); % Define Gamma(-1), i.e. the cost, for each state
%cost = ones(Ns, 1) * Inf;

%statemap = (1:Ns).';

%cool but useless:
%newStateBase=@(currstate) (mod(currstate-1, M^(N1+N2-1)) * 4);


% --- Main loop

for k = 1 : length(r)
    
    % "Pointer" to the column we are gonna write in this iteration
    survSeq_writingcol = survSeq_writingcol + 1;
    
    % Initialize the costs of the new states to -1
    costnew = - ones(Ns, 1);
    
    % Vector of the predecessors: the i-th element is the predecessor at
    % time k-1 of the i-th state at time k.
    pred = zeros(Ns, 1);
    
    % Counter for the new state. It is determined iteratively even though a
    % closed form expression exists: (mod(state-1, M^(N1+N2-1)) * M + j).
    newstate = 0;
    
    for state = 1 : Ns  % Cycle through all states, at time k-1 (?)
        
        for j = 1:M     % M possibilities for the new symbol
            
            % Index of the new state: it's mod(state-1, M^(N1+N2-1)) * M + j
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            
            % TODO optimize
            % Supposed new sequence, obtained from the old state adding a new symbol
            supposednewseq = [symb(mod(survSeq(state, 1:survSeq_writingcol-1)-1,M)+1), symb(j)];
            %supposednewseq = [symb(mod(survSeq(statemap(state), 1:survSeq_writingcol-1)-1,M)+1), symb(j)];
            
            % TODO optimize
            % Compute desired signal u assuming the input sequence is the one above
            difflength = N - length(supposednewseq);
            if difflength > 0
                supposednewseq = [zeros(1, difflength), supposednewseq];
            else
                supposednewseq = supposednewseq(end-N+1:end);
            end
            u = supposednewseq * flipud(hi);
            
            % Compute the cost of the new state assuming this input sequence,
            % then update the cost of the new state, and overwrite the predecessor,
            % if this transition has a lower cost than before.
            newstate_cost = cost(state) + abs(r(k) - u)^2;
            if costnew(newstate) == -1 ...     % not assigned yet, or...
                    || costnew(newstate) > newstate_cost    % ...found path with lower cost
                % Update everything
                costnew(newstate) = newstate_cost;
                pred(newstate) = state;
            end
            
            % fprintf('k=%d\t%d->%d\tu=%f+j%f\tr=%f+j%f\n', k, state, newstate, real(u), imag(u), real(r(k)), imag(r(k)))
            
        end
    end
    
    % Handle memory. If we were gonna write to the last column, first we store
    % the oldest chunk of the decided sequence (we decide it for good), then
    % we shift the matrix to make room for new data.
    if survSeq_writingcol == size(survSeq, 2)
        
        % Take the first row of old decided states and store it before erasing it.
        detectedStates(1+survSeq_shift : survSeq_shift+MEMORY) = survSeq(1, 1:MEMORY);
        
        % Shift of half the size, that is our memory.
        survSeq(:, 1:MEMORY) = survSeq(:, MEMORY+1:end);
        survSeq_shift = survSeq_shift + MEMORY;
        survSeq_writingcol = MEMORY; % Actually it is (SIZE-MEMORY)
    end
    
        
    % The following operations strongly affect the computation time, if
    % the number of states is not too large. Otherwise the bottleneck is
    % the number of iterations of the inner loops above.
    % TODO optimize the way matrices are handled. Avoid creating so many.
    
    temp = zeros(size(survSeq));
    for newstate = 1:Ns
        temp(newstate, 1:survSeq_writingcol) = ...
            [survSeq(pred(newstate), 1:survSeq_writingcol-1), newstate];
    end
    survSeq = temp;
    
    
%     % For each new state, see where its predecessor is (which row according to
%     % statemap), then add newstate at the end of that survival sequence, and iterate
%     % for all new states. Finally, update the statemap according to the new states at k.
%     % Instead of sorting the rows according to the new states at k, we
%     % store a map so as to keep track of the row in which a certain state
%     % at time k is.
%     
%     % These are the states at k-1 that were discarded by the algorithm, so
%     % we can overwrite them
%     deadsequences = setdiff(1:Ns, unique(pred)); % data in A that is not in B
%     deadseq_idx = 1;
%     alreadyused = false(Ns, 1); % refers to the true index of the matrix
%     for newstate = 1:Ns
%         if ~alreadyused(statemap(pred(newstate)))
%             survSeq(statemap(pred(newstate)), survSeq_writingcol) = newstate;
%             alreadyused(statemap(pred(newstate))) = true;
%         else
%             survSeq(statemap(deadsequences(deadseq_idx)), 1:survSeq_writingcol) = ...
%                 [survSeq(statemap(pred(newstate)), 1:survSeq_writingcol-1), newstate];
%             alreadyused(statemap(deadsequences(deadseq_idx))) = true; % shouldn't need this, I think
%             deadseq_idx = deadseq_idx + 1;
%         end
%     end
%     statemap = survSeq(:, survSeq_writingcol);
    
    
    
    %allcosts(:, k) = cost;
    
    % Update the cost to be used as cost at time k-1 in the next iteration
    cost = costnew;
end

toc
elapsed_time = toc;



% --- Finish storing, then get the symbols

detectedStates(1+survSeq_shift : survSeq_shift+survSeq_writingcol-1) = ...
    survSeq(1, 1:survSeq_writingcol-1);
detected = symb(mod(detectedStates-1, M) + 1);



% --- Output results

num_errors = sum(packet-detected.' ~= 0);
[pbit, num_bit_errors] = BER(packet, detected);

fprintf('P_err is approximately %.g (%d errors)\n', num_errors / length(packet), num_errors)

% figure, hold on
% stem(real(r))
% stem(real(packet), 'x')
% stem(real(detected), 'd')
% legend('r', 'sent', 'detected')
% figure, stem(0:N1+N2, real(hi)), hold on, stem(0:N1+N2, imag(hi)), title('hi')

end