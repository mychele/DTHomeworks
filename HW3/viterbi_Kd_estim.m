function [iecdfyavg, ecdfx] = viterbi_Kd_estim()

% Runs numsim simulations of the transmission of 3 ML sequences with length 2^20 - 1, with
% a Viterbi detector that has Kd = MAX_Kd. At each iteration (starting after the transient
% Kd), it checks which is the minimum value of Kd that would yield no premature decision.
% That is, all survivor sequences contain the same symbols from a certain column
% backwards, and that column is given by the minimum Kd needed, which is a sort of
% "memory" of the Viterbi algorithm.
% If the resulting complementary CDF is y for a certain value of Kd, it means that for that
% Kd, the number of iterations in which the decision on the oldest symbol is actually made
% among different symbols is y times the total number of iterations (except for the
% transient in which the survivor sequence is filled for the first time). In those cases,
% the symbol that belongs to the most likely sequence is chosen, which might turn out to
% be a suboptimal decision. The number of times we make this approximation on the Viterbi
% algorithm is y times the length of the sequence.

L_data = 3*(2^20 - 1);

numsim = 16;
MAX_Kd = 80;
snr_channel = 5;   % in dB, worst case scenario in our simulations
ecdfy = zeros(numsim, MAX_Kd+1);

% From exercise 1
N1 = 0;
N2 = 4;
assumed_dly = 2;
assumed_m_opt = 10;

% Peform numsim simulation using 16 parallel threads
parpool(16);
parfor sim_i = 1:numsim
    
    % --- Create, send and receive data, estimate channel and prepare for detection
    
    % Create, send and receive data with the given channel
    [packet, r, ~] = txrc(L_data, snr_channel, assumed_m_opt);
    
    % Estimate the channel using the first 100 samples (4*length(ts))
    [ h, ~ ] = get_channel_info(r(assumed_dly+1:25+assumed_dly), N1, N2);
    
    % Normalize x and h
    x = r / h(N1+1).';    % data normalized by h0
    hi = h / h(N1+1).';   % impulse response normalized by h0
    
    
    % --- Perform detection and get bin counts
    [ecdfy_this] = viterbi_Kd_estim_helper( ...
        packet, x(1+assumed_dly-N1:end), hi, N1, N2, 0, N2, MAX_Kd);
    ecdfy(sim_i, :) = ecdfy_this;
end


delete(gcp);  % Delete parallel pool


% Get and store results
iecdfyavg = 1 - mean(ecdfy);
ecdfx = 0:MAX_Kd;
save('viterbi_Kd_estimation', 'iecdfyavg', 'ecdfx', 'numsim');

% Plot
% stairs(ecdfx, iecdfyavg)
% title(['Complementary CDF averaged over ', int2str(numsim), ' simulations'])
% xlabel('K_{d, min}')
% set(gca, 'YScale', 'log')
% ylim([1e-8 1e-3]), xlim([find(iecdfyavg < 1e-3, 1)-1, find(iecdfyavg < 1e-8, 1)+1])


end



function [ ecdfy ] = viterbi_Kd_estim_helper( packet, r, hi, N1, N2, L1, L2, MAX_Kd )
% Returns the values of the empirical CDF of Kd_min, for Kd_min = 0, ..., MAX_Kd, given
% the usual Viterbi parameters.

if (L1 > N1) || (L2 > N2)
    disp('The considered precursors and postcursors cannot be more that the actual ones!')
    return
end


% --- Setup

M = 4;
symb = [1+1i, 1-1i, -1+1i, -1-1i]; % All possible transmitted symbols (QPSK)
Ns = M ^ (L1+L2); % Number of states
r  =  r(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of r
hi = hi(1+N1-L1 : end-N2+L2);   % Discard initial and final samples of hi



% ------------------
% --- Init stuff ---
% ------------------

tStart = tic;   % Use a variable to avoid conflict with parallel calls to tic/toc
survSeq = zeros(Ns, MAX_Kd);
detectedSymb = zeros(1, length(packet));
cost = zeros(Ns, 1); % Define Gamma(-1), i.e. the cost, for each state

Kd_min = zeros(1, length(packet)-MAX_Kd);


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
        temp(newstate, 1:MAX_Kd) = ...    % In the new matrix except the last col
            [survSeq(pred(newstate), 2:MAX_Kd), ... % we write the data we had except the oldest one,
            symb(mod(newstate-1, M)+1)];        % and then the new symbol we just supposed to have received.
    end
    [~, decided_index] = min(costnew);      % Find the oldest symbol that yields the min cost
    detectedSymb(1+k) = survSeq(decided_index, 1); % and store it (decide it for good).
    survSeq = temp;
    
    % Update the cost to be used as cost at time k-1 in the next iteration
    cost = costnew;
    
    % ------------------- This is specific to the Kd estimation function
    %  Check what would be, at this iteration k, the minimum Kd in order to not discard
    %  any useful information. Store this value, compute the ecdf at the end.
    %  This is only done from k==Kd+1 on.
    if k > MAX_Kd
        for i = 0 : MAX_Kd-1
            if length(unique(survSeq(:, MAX_Kd-i))) == 1
                Kd_min(k-MAX_Kd) = i;
                break;
            end
        end
    end
end

toc(tStart)



% --- Output results

[ecdfy, ecdfx] = ecdf(Kd_min);
bincenters = 0:MAX_Kd;
[bincount, ~] = ecdfhist(ecdfy, ecdfx, bincenters);
ecdfy = cumsum(bincount);

end