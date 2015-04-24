% Impulse response estimation
% We are at the receiver, we know what the sender is sending and we try to
% estimate it with the LS method (for reference, see page 244).

% Note: As the receiver, we do _not_ know neither N_h nor sigma_w

clear, clc, close all
rng default

%% Generate gi
channel_generator;

%% Loop to determine suitable values of N, L

printmsg_delete = '';
maxN = 10;

% time counter that let the desired output d to be computed with a
% different impulse response at every iteration, as it would happen in
% reality. 
time = 1;

for L = [3, 7, 15, 31]
    for N = 1:maxN % Supposed length of the impulse response of the channel
        printmsg = sprintf('L = %d, N = %d\n', L, N);
        fprintf([printmsg_delete, printmsg]);
        printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
        
        % Supposed length of the impulse response of the channel in each polyphase branch.
        n_short = mod(4-N, 4); % Num branches with a shorter filter than others
        % N_i is the number of coefficients of the filter of the i-th branch.
        N_i(1:4-n_short) = ceil(N/4);
        N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
        
        %% Generate training sequence
        % The x sequence must be a partially repeated M-L sequence of length L. We
        % need it to have size L+N-1.
        % To observe L samples, we need to send L+N-1 samples of the training
        % sequence {x(0), ..., x((N-1)+(L-1))}
        p = MLsequence(L);
        x = [p; p(1:max(N_i)-1)];
        x(x == 0) = -1;
        
        %% Estimation of h and d multiple times
        numsim = 300; % It seems to converge even with small values of numsim,
        % lowered down to 200 in order to make computation feasible with a
        % veery big g_mat
        error_func_temp = zeros(numsim, 1);
        for k =1:numsim
            [d, h_mean] = channel_output(x, T, Tc, sigma_w, N_h, g_mat(:, time:end));
            time = time + (L+N)*T/Tc; %(L+N)*4, they shouldn't overlap and 
            % there should be enough impulse responses. Probably we need
            % less!
            [h_hat, d_hat] = h_estimation(x, d, L, N_i);
            d_no_trans = d(end-length(d_hat)+1 : end);
            error_func_temp(k) = sum(abs(d_hat - d_no_trans).^2);
        end
        error_func(N) = mean(error_func_temp);
    end
    
    % NOTE: the receiver doesn't know the reference value to which the
    % functional should tent to, it is plotted for debugging purposes
    figure, plot(10*log10(error_func)), hold on, plot(1:N, ones(1, N)*10*log10(sigma_w*length(d_hat)))
    xlabel('N')
    ylabel('\epsilon [dB]')
    legend('Error functional Eps', 'Theoretical value of Eps')
    grid on, title(['Error function with L=', int2str(L)])
end

return

%% Plot h and h_hat

figure, stem(abs(h_mean), 'DisplayName', 'First coefficient of each poly branch (mean)')
legend('-DynamicLegend'), hold on,
for k = 2:N_i
    stem(zeros(4,1), 'DisplayName', strcat(num2str(k), ' coefficient of each poly branch'))
    legend('-DynamicLegend'), hold on,
end
for k = 1:N_i
    stem(abs(h_hat(:, k)), 'DisplayName', 'hhat')
    legend('-DynamicLegend')
end
ax = gca; ax.XTick = 1:4;
xlim([0.5, 4.5])

%% Compute estimation error |h - h_hat|^2

% just to ease computation set to 0 the unused coefficients
% reshape hhat and h_mean
h_hat_array = reshape(h_hat, 4*max(N_i), 1); % this vector represents
% the estimated hhat_i for n = 0, 1, ... N - 1. Its length is actually more
% than N, but the last elements are just 0, therefore they can be removed
h_hat_array = h_hat_array(1:N);
h_mean_array = h_mean(1:N_h);
% the vector to which I compare the estimate has length N_h, I should
% resize the shortest vector to have the same length of the other by adding
% some 0
if N < N_h
    h_hat_array = [h_hat_array; zeros(N_h - N, 1)];
elseif N > N_h
    h_mean_array = [h_mean_array; zeros(N - N_h, 1)];
end % if N=N_h already ok
errorpower = sum(abs(h_hat_array - h_mean_array).^2);



%% Repeat the estimate 1000 times (note, this will be influenced by the actual
% choice of N and L
numsim = 1000;

hhat_mat = zeros(4*max(N_i), numsim);
errorpower_array = zeros(1, numsim);
for simiter = 1:numsim
    disp(simiter);
    
    % it would be better to generate d offline?
    d = zeros(T*length(x),1);
    g_used_coeff = zeros(4, length(x)); % check these dimensions
    for k = 0:length(x)-1
        % Generate white noise
        w = wgn(4, 1, 10*log10(sigma_w), 'complex');
        for idx = 0:4-1 % actually idx varies between 0 and 4 since there are four branches
            if (idx < N_h)
                d(k*T + idx*Tc + 1) = g_mat(idx+1,k*T + idx*Tc+1) * x(k+1) + w(idx+1);
            else
                d(k*T + idx*Tc + 1) = 0 + w(idx+1); % no ray, just the noise
            end
        end
    end
    
    % create four different d_i vector, by sampling at step 4 the complete vector
    % d. Each of them is the output of the polyphase brach at "lag" i
    d_poly = zeros(length(d)/4, 4); % each column is a d_i
    for idx  = 1:4
        d_poly(:, idx) = d(idx:4:end);
    end
    
    % Using the data matrix (page 246), easier implementation
    h_hat = zeros(4,max(N_i)); % estimate 4 polyphase represantations
    % this matrix is dimensioned to have the maximum number of coefficients
    % for each of the four branch, the unused (i.e. unestimated) ones will be
    % left zero
    for idx = 1:4
        if N_i(idx) > 0
            I = zeros(L,N_i(idx));
            for column = 1:N_i(idx)
                I(:,column) = x(N_i(idx)-column+1:(N_i(idx)+L-column));
            end
            o = d_poly(N_i(idx):N_i(idx) + L - 1, idx);
            
            Phi = I'*I;
            theta = I'*o;
            
            h_hat(idx, 1:N_i(idx)) = inv(Phi) * theta;
        end % if N_branch is 0 don't estimate and leave hhat to 0
    end
    hhat_mat(:, simiter) = reshape(h_hat, 4*max(N_i), 1);
    
    
    % just to ease computation set to 0 the unused coefficients
    % reshape hhat and h_mean
    h_hat_array = reshape(h_hat, 4*max(N_i), 1); % this vector represents
    % the estimated hhat_i for n = 0, 1, ... N - 1. Its length is actually more
    % than N, but the last elements are just 0, therefore they can be removed
    h_hat_array = h_hat_array(1:N);
    h_mean_array = h_mean(1:N_h);
    % the vector to which I compare the estimate has length N_h, I should
    % resize the shortest vector to have the same length of the other by adding
    % some 0
    if N < N_h
        h_hat_array = [h_hat_array; zeros(N_h - N, 1)];
    elseif N > N_h
        h_mean_array = [h_mean_array; zeros(N - N_h, 1)];
    end % if N=N_h already ok
    errorpower_array(simiter) = sum(abs(h_hat_array - h_mean_array).^2);
end

hhat_array_est = mean(hhat_mat, 2);
% reshape it
hhat_mat_mean = reshape(hhat_array_est, 4, max(N_i));


figure, stem(abs(h_mean), 'DisplayName', 'First coefficient of each poly branch (mean)')
legend('-DynamicLegend'), hold on,
for k = 2:N_i
    stem(zeros(4,1), 'DisplayName', strcat(num2str(k), ' coefficient of each poly branch'))
    legend('-DynamicLegend'), hold on,
end
for k = 1:N_i
    stem(abs(hhat_mat_mean(:, k)), 'DisplayName', 'hhat')
    legend('-DynamicLegend')
end
ax = gca; ax.XTick = 1:4;
xlim([0.5, 4.5])
title('mean of hhat over 1000 estimates')

% an estimate of E(|hhat - h|^2) over 1000 realizations could be given by
% the mean of E at each iteration
errorpower_est = mean(errorpower_array);
