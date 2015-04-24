%% Impulse response estimation
% We are at the receiver, we know what the sender is sending and we try to
% estimate it with the LS method (for reference, see page 244).

% Note: As the receiver, we do _not_ know neither N_h nor sigma_w

clear all
close all
clc

%% Generate gi

N_h = 3;
a_dopp = [1, -4.4153, 8.6283, -9.4592, 6.1051, -1.3542, -3.3622, 7.2390, ...
    -7.9361, 5.1221, -1.8401, 2.8706e-1];
b_dopp = [1.3651e-4, 8.1905e-4, 2.0476e-3, 2.7302e-3, 2.0476e-3, 9.0939e-4, ...
    6.7852e-4, 1.3550e-3, 1.8076e-3, 1.3550e-3, 5.3726e-4, 6.1818e-5, -7.1294e-5, ...
    -9.5058e-5, -7.1294e-5, -2.5505e-5, 1.3321e-5, 4.5186e-5, 6.0248e-5, 4.5186e-5, ...
    1.8074e-5, 3.0124e-6];
[h_dopp, ~] = impz(b_dopp, a_dopp);
hds_nrg = sum(h_dopp.^2);
b_dopp = b_dopp / sqrt(hds_nrg);
Tc = 1; % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Kdb = 3; % 3 dB
K = 10^(Kdb/10);
fd = 5*10^-3/T; % doppler frequency
Tq = Tc; % Fundamental sampling time. This is the same as Tc, right???
% We stick to anastasopolous chugg paper (1997) and choose Tp such that
% fd*Tp = 0.1
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
tau = 0:Tc:N_h-1;
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise
pdp_gauss = 1/tau_rms * exp(-tau/tau_rms);

C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|gtilde_i|^2]) = 1 - C^2
pdp_gauss = pdp_gauss.*(1-C^2)/sum(pdp_gauss);

M_d = sum(pdp_gauss);

snr = 10; % db
snr_lin = 10^(snr/10);
sigma_w = 1/(T/Tc*snr); % the PN sequence has power 1


fprintf('fd * Tp = %d \n', Tp*fd);
g_samples_needed = 200000; % Some will be dropped because of transient
w_samples_needed = ceil(g_samples_needed / Tp);
transient = ceil(g_samples_needed/4);

% Generate complex-valued Gaussian white noise with zero mean and unit
% variance
%rng('default');
g_mat = zeros(N_h, g_samples_needed - transient);
for ray = 1:N_h
    w = wgn(w_samples_needed,1,0,'complex');
    fprintf('variance of white noise=%d \n', var(w))
    % Filter the wgn with a narrowband filter. The filter will have the
    % classical Doppler spectrum in the frequency domain, with f_d * Tc = 1.25*10^-3
    % By using the approach suggested in anastchugg97 we use an iir filter
    % with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
    % Note that it is
    % possible to initialize properly the filters, consider doing it in
    % following versions since it allows to start in steady state conditions
    % and avoid dropping about many samples (it is a very cool thing!)
    gprime = filter(b_dopp, a_dopp, w);
    fprintf('g %d after iterpolation has mean %d and variance %d  \n', ray, mean(gprime), var(gprime));
    
    % Interpolation
    t = 1:length(gprime);
    t_fine = Tq/Tp:Tq/Tp:length(gprime);
    
    g_fine = interp1(t, gprime, t_fine, 'spline');
    fprintf('g %d after iterpolation has mean %d and variance %d  \n', ray, mean(g_fine), var(g_fine));
    % Drop the transient
    g_mat(ray, :) = g_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h
    g_mat(k, :) = g_mat(k, :)*sqrt(pdp_gauss(k));
end

% Only for LOS component
g_mat(1, :) = g_mat(1, :) + C;

clear a_dopp b_dopp g_fine gprime h_dopp pdp_gauss t_dopp t_fine g_samples_needed w_samples_needed

%% Actual simulation

% To observe L samples, we need to send L+N-1 samples of the training
% sequence {x(0), ..., x((N-1)+(L-1))}
% TODO Try multiple combinations. Generally, L = 2*N
L = 15; % Length of the observation
% Supposed length of the impulse response of the channel
N_max = 12;
errorpower_est = zeros(1, N_max);
for N = 1:N_max % Supposed length of the impulse response of the channel
    
    fprintf('Analyzing the channel, supposed length %d \n', N)
    % Supposed length of the impulse response of the channel in each polyphase branch
    % the number of branches that have a coefficient less than the other
    % branches is
    n_short = mod(4-N, 4);
    % then the number of coefficients in each branch is
    N_branch(1:4-n_short) = ceil(N/4);
    N_branch(4-n_short + 1:4) = ceil(N/4) - 1;
    
    
    %% Generate training sequence
    % The x sequence must be a partially repeated M-L sequence of length L. We
    % need it to have size L+N-1.
    % this holds for L = 15
    r = log2(L+1);
    p = zeros(L,1);
    p(1:r) = ones(1,r).'; % Set arbitrary initial condition
    for l = r+1:(L)
        p(l) = xor(p(l-3), p(l-4)); % -1, -2 for L=3
    end
    clear l
    
    x = [p; p(1:N_branch-1)];
    x(x == 0) = -1;
    
    %% Generate a first output in order to track which coefficients are actually used
    d = zeros(T*length(x),1);
    g_used_coeff = zeros(4, length(x)); % check these dimensions
    for k = 0:length(x)-1
        % Generate white noise
        w = wgn(4, 1, 10*log10(sigma_w), 'complex');
        for idx = 0:4-1 % actually idx varies between 0 and 4 since there are four branches
            if (idx < N_h)
                d(k*T + idx*Tc + 1) = g_mat(idx+1,k*T + idx*Tc+1) * x(k+1) + w(idx+1);
                % store the coefficient actually used, it will be useful later on
                g_used_coeff(idx + 1, k + 1) = g_mat(idx+1,k*T + idx*Tc+1);
            else
                d(k*T + idx*Tc + 1) = 0 + w(idx+1); % no ray, just the noise
                g_used_coeff(idx + 1, k + 1) = 0;
            end
        end
    end
    
    % remove the coefficients of the transient from the g_used_coeff matrix
    g_used_coeff = g_used_coeff(:, max(N_branch) - 1 + 1:end); % N-1 is the max transient, +1 because of MATLAB
    
    h_mean = mean(g_used_coeff, 2);
    
    %% Repeat the estimate 1000 times (note, this will be influenced by the actual
    % choice of N and L
    numsim = 1000;
    
    hhat_mat = zeros(4*max(N_branch), numsim);
    errorpower_array = zeros(1, numsim);
    for simiter = 1:numsim
        %fprintf(num2str(simiter));
        
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
        h_hat = zeros(4,max(N_branch)); % estimate 4 polyphase represantations
        % this matrix is dimensioned to have the maximum number of coefficients
        % for each of the four branch, the unused (i.e. unestimated) ones will be
        % left zero
        for idx = 1:4
            if N_branch(idx) > 0
                I = zeros(L,N_branch(idx));
                for column = 1:N_branch(idx)
                    I(:,column) = x(N_branch(idx)-column+1:(N_branch(idx)+L-column));
                end
                o = d_poly(N_branch(idx):N_branch(idx) + L - 1, idx);
                
                Phi = I'*I;
                theta = I'*o;
                
                h_hat(idx, 1:N_branch(idx)) = inv(Phi) * theta;
            end % if N_branch is 0 don't estimate and leave hhat to 0
        end
        hhat_mat(:, simiter) = reshape(h_hat, 4*max(N_branch), 1);
        
        
        % just to ease computation set to 0 the unused coefficients
        % reshape hhat and h_mean
        h_hat_array = reshape(h_hat, 4*max(N_branch), 1); % this vector represents
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
    hhat_mat_mean = reshape(hhat_array_est, 4, max(N_branch));
    
    % an estimate of E(|hhat - h|^2) over 1000 realizations could be given by
    % the mean of E at each iteration
    errorpower_est(N) = mean(errorpower_array);
    fprintf('\n');
end

figure, plot(errorpower_est)
xlabel('N which tracks the real N_h')
ylabel('Estimate of E(|h - hhat|^2) across 1000 realizations')
ax = gca; ax.XTick = 1:N_max;

figure, plot(10*log10(errorpower_est))
xlabel('N which tracks the real N_h')
ylabel('Estimate of E(|h - hhat|^2) across 1000 realizations [dB]')
ax = gca; ax.XTick = 1:N_max;
