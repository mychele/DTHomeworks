function [ detected, pbit, num_bit_errors ] = viterbi( packet, r, hi, N1, N2, L1, L2, traininglength )
%VITERBI

% State: most recent is at the left, as in the book, and it has lowest weight
% (thus, in the base-M representation, it is of course the rightmost one).

%fprintf('Viterbi started.\n')
if (L1 > N1) || (L2 > N2)
    disp('The considered precursors and postcursors cannot be more that the actual ones!')
    return
end


% --- Setup

M = 4;
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % All possible transmitted symbols (QPSK)
N = N1 + N2 + 1;
Kd = 20 * N;      % Trellis size, that is the size of the state matrix as well
Ns = M ^ (L1+L2); % Number of states
r  =  r(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of r
hi = hi(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of hi



% ------------------
% --- Init stuff ---
% ------------------

tStart = tic;   % Use a variable to avoid conflict with parallel calls to tic/toc
survSeq = zeros(Ns, Kd);
%survSeq(:, Kd+1) = symb(mod(0:Ns-1, M) + 1);
% TODO: init the first elements to the end of ML sequence.
detectedSymb = zeros(1, length(packet));
cost = zeros(Ns, 1); % Define Gamma(-1), i.e. the cost, for each state


% -- Define u_mat matrix

statelength = L1 + L2; % number of digits of the state, i.e. its length in base M=4
statevec = zeros(1, statelength); % symbol index, from the oldest to the newest
u_mat = zeros(Ns, M);
for state = 1:Ns
    
    % Set value of the current element of u_mat
    for j = 1:M
        lastsymbols = [symb(statevec + 1), symb(j)]; % symbols, from the oldest to the newest
        u_mat(state, j) = lastsymbols * flipud(hi);
    end
    
    % Update statevec
    statevec(statelength) = statevec(statelength) + 1;
    i = statelength;
    while (statevec(i) >= M && i > 1)
        statevec(i) = 0;
        i = i-1;
        statevec(i) = statevec(i) + 1;
    end
end




% -----------------
% --- Main loop ---
% -----------------

for k = 1 : length(r)
    
    % Initialize the costs of the new states to -1
    costnew = - ones(Ns, 1);
    
    % Vector of the predecessors: the i-th element is the predecessor at
    % time k-1 of the i-th state at time k.
    pred = zeros(Ns, 1);
    
    % Counter for the new state. It is determined iteratively even though a
    % closed form expression exists: (mod(state-1, M^(N1+N2-1)) * M + j).
    newstate = 0;
    
    
    for state = 1 : Ns  % Cycle through all states, at time k-1
        
        for j = 1 : M   % M possibilities for the new symbol
            
            % Index of the new state: it's mod(state-1, M^(N1+N2-1)) * M + j
            newstate = newstate + 1;
            if newstate > Ns, newstate = 1; end
            
            % Desired signal u assuming the input sequence is the one given by the current
            % state "state", followed by a new assumed symbol given by j.
            u = u_mat(state, j);
            
            % Compute the cost of the new state assuming this input sequence, then update
            % the cost of the new state, and overwrite the predecessor, if this transition
            % has a lower cost than before.
            newstate_cost = cost(state) + abs(r(k) - u)^2;
            if costnew(newstate) == -1 ...     % not assigned yet, or...
                    || costnew(newstate) > newstate_cost  % ...found path with lower cost
                costnew(newstate) = newstate_cost;
                pred(newstate) = state;
            end
        end
    end
    
    
    % Update the survivor sequence by shifting the time horizon of the matrix by one, and
    % rewrite the matrix with the new survival sequences sorted by current state.
    % Meanwhile, decide the oldest sample (based on minimum cost) and get rid of it to
    % keep only Kd columns in the matrix.
    temp = zeros(size(survSeq));
    for newstate = 1:Ns
        temp(newstate, 1:Kd) = ...    % In the new matrix except the last col
            [survSeq(pred(newstate), 2:Kd), ... % we write the data we had except the oldest one,
            symb(mod(newstate-1, M)+1)];        % and then the new symbol we just supposed to have received.
    end
    [~, decided_index] = min(costnew);      % Find the oldest symbol that yields the min cost
    detectedSymb(1+k) = survSeq(decided_index, 1); % and store it (decide it for good).
    survSeq = temp;
    
    % Update the cost to be used as cost at time k-1 in the next iteration
    cost = costnew;
end

toc(tStart)



% --- Finish storing, then get the symbols
detectedSymb(length(r)+2 : length(r)+Kd) = survSeq(decided_index, 1:Kd-1);
% Decided using the min cost from the last iteration
detectedSymb = detectedSymb(Kd+1 : end);
detected = detectedSymb;
detected = detected(2:length(packet)+1);    % Discard first symbol (time k=-1)
detected = detected (1+traininglength : end);  % Discard detected training sequence both
packet   = packet   (1+traininglength : end);  % in the received and in the detected data



% --- Output results

[pbit, num_bit_errors] = BER(packet, detected);

num_errors = sum(packet-detected.' ~= 0);
fprintf('P_err = %.g (%d errors)\n', num_errors / length(packet), num_errors)
fprintf('P_bit = %.g (%d errors)\n', pbit, num_bit_errors)

% > Very useful plot for debugging <
% figure, hold on
% stem(-L1-traininglength:length(r)-1-L1-traininglength, real(r))
% stem(0:length(packet)-1, real(packet), 'x')
% stem(0:length(detected)-1, real(detected), 'd')
% legend('r', 'sent', 'detected')

end