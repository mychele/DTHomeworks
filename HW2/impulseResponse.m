% Impulse response estimation
% We are at the receiver, we know what the sender is sending and we try to
% estimate it with the LS method (for reference, see page 244).
% Note: As the receiver, we do _not_ know neither N_h nor sigma_w

clear, clc, close all
rng default

%% Generate time-variant i.r. of the channel and initialize everything
channel_generator;
sigma_w = 1/(T/Tc*snr); % the PN sequence has power 1

%% Loop to determine suitable values of N, L

printmsg_delete = ''; % Just to display progress updates
maxN = 10;
% Note that with maxN<13 we don't have problems with the condition N<=L.

% Time counter that allows the output d to be computed with a different
% impulse response at every iteration, as it would happen in reality.
time = 1;
L_vec = [3, 7, 15, 31, 63, 127];
numsim = 100; % It seems to converge even with small values of numsim
error_func = zeros(length(L_vec), maxN, numsim);
for L_index = 1:length(L_vec)
    L = L_vec(L_index);
    
    % --- Generate training sequence
    % The x sequence must be a partially repeated M-L sequence of length L. We
    % need it to have size L+N-1. To observe L samples, we need to send L+N-1
    % samples of the training sequence {x(0), ..., x((N-1)+(L-1))}.
    p = MLsequence(L);
    x = [p; p(1:ceil(maxN/4)-1)]; % create a seq which is long enough for the maximum N
    x(x == 0) = -1;
    
    % --- Estimation of h and d multiple times
    for k =1:numsim
        
        % Print progress update
        printmsg = sprintf('L = %d, simulation number %d\n', L, k);
        fprintf([printmsg_delete, printmsg]);
        printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
        
        % Transmit only one time and estimate h for different N
        [d, ~] = channel_output(x, T, Tc, sigma_w, N_h, h_mat(:, time:end));
        time = time + 50*(L+maxN)*T/Tc; % the time windows are sufficiently spaced apart
        for N = 1:maxN % N is the supposed length of the impulse response of the channel
            % Compute the supposed length of each branch
            n_short = mod(4-N, 4); % Num branches with a shorter filter than others
            % N_i is the number of coefficients of the filter of the i-th branch.
            N_i(1:4-n_short) = ceil(N/4);
            N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
            [h_hat, d_hat] = h_estimation( x(end-(L+max(N_i)-1) + 1 : end), ...
                        d(end - 4*(L+max(N_i)-1) + 1: end), L, N_i);
            d_no_trans = d(end-length(d_hat)+1 : end);
            error_func(L_index, N, k) = sum(abs(d_hat - d_no_trans).^2)/length(d_hat);
        end
    end
end
error_func = mean(error_func, 3);

% Plot the empirical error functional for different pairs (L, N)
figure, hold on
for i = 1:length(L_vec)
    plot(10*log10(error_func(i, :)), 'DisplayName', strcat('L=', num2str(L_vec(i))))
    legend('-DynamicLegend')
end
xlabel('N'), ylabel('\epsilon [dB]'), title('Error function')
grid on, box on, ylim([-20, -10])

%% Estimate E(|h-hhat|^2) by repeating the estimate 1000 times and assuming
% h known

printmsg_delete = '';

% time counter that let the desired output d to be computed with a
% different impulse response at every iteration, as it would happen in
% reality.
time = 1;
L_vec = [3, 7, 15, 31];
numsim = 1000;
deltah_square = zeros(length(L_vec), numsim, maxN);
for L_index = 1:length(L_vec)
    L = L_vec(L_index);
    
    % --- Generate training sequence
    % The x sequence must be a partially repeated M-L sequence of length L. We
    % need it to have size L+N-1. To observe L samples, we need to send L+N-1
    % samples of the training sequence {x(0), ..., x((N-1)+(L-1))}.
    p = MLsequence(L);
    x = [p; p(1:ceil(maxN/4)-1)]; % create a seq which is long enough for the maximum N
    x(x == 0) = -1;
    
    % --- Estimation of h multiple times
    for k=1:numsim
        printmsg = sprintf('L = %d, simulation number = %d\n', L, k);
        fprintf([printmsg_delete, printmsg]);
        printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
        
        % We transmit only one time and then estimate h for different N
        [d, h_mean] = channel_output(x, T, Tc, sigma_w, N_h, h_mat(:, time:end));
        time = time + 50*(L+maxN)*T/Tc; % the time windows are sufficiently spaced apart
        
        for N = 1:maxN % N is the supposed length of the impulse response of the channel
            n_short = mod(4-N, 4); % Num branches with a shorter filter than others
            % N_i is the number of coefficients of the filter of the i-th branch.
            N_i(1:4-n_short) = ceil(N/4);
            N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
            
            % LS estimation of h
            [h_hat, ~] = h_estimation(x, d, L, N_i);
            
            % Compute delta_h squared
            h_hat_array = reshape(h_hat, 4*max(N_i), 1);
            h_hat_array = h_hat_array(1:N);
            h_mean_array = h_mean(1:N_h);
            % The vector to which we compare the estimate has length N_h,
            % the estimated h_hat has length N. Now we make them the same length.
            if N < N_h
                h_hat_array = [h_hat_array; zeros(N_h - N, 1)];
            elseif N > N_h
                h_mean_array = [h_mean_array; zeros(N - N_h, 1)];
            end % if N=N_h already ok
            deltah_square(L_index, N, k) = sum(abs(h_hat_array - h_mean_array).^2);
        end
    end
end
deltah_square = mean(deltah_square, 3);

%% Compare with theoretical

deltah_square_theor = zeros(length(L_vec), maxN);
for L_index = 1:length(L_vec)
    L = L_vec(L_index);
    for N = 1:maxN
        if(ceil(N/4) > L)   % The estimate cannot be performed (see report)
            deltah_square_theor(L_index, N) = NaN;
            break
        end
        
        n_short = mod(4-N, 4); % Num branches with a shorter filter than others
        % N_i is the number of coefficients of the filter of the i-th branch.
        N_i(1:4-n_short) = ceil(N/4);
        N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
        
        deltah_square_theor(L_index, N) = ...
                sigma_w / (L+1) * sum(N_i .* (L+2-N_i) ./ (L+1-N_i));
    end
end

% Plot results
figure, hold on
for L_index = 1:4
    plot(10*log10(deltah_square(L_index, :)), ...
            'DisplayName', sprintf('L=%d experimental', L_vec(L_index)))
    plot(10*log10(deltah_square_theor(L_index, :)),'-.', ...
            'DisplayName', sprintf('L=%d theoretical', L_vec(L_index)))
    legend('-DynamicLegend')
end
xlabel('N that tracks the real N_h'), ylabel('Estimate of E(|h - hhat|^2) [dB]')
title('Estimate of E(|h - hhat|^2) across 1000 realizations')
ax = gca; ax.XTick = 1:maxN;
ylim([-30 -5]), grid on, box on