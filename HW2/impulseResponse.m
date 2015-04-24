%% Impulse response estimation
% We are at the receiver, we know what the sender is sending and we try to
% estimate it with the LS method (for reference, see page 244).

% Note: As the receiver, we do _not_ know neither N_h nor sigma_w

clear
close all
clc

% To observe L samples, we need to send L+N-1 samples of the training
% sequence {x(0), ..., x((N-1)+(L-1))}
% TODO Try multiple combinations. Generally, L = 2*N
L = 15; % Length of the observation
N = 10; % Supposed length of the impulse response of the channel
% Supposed length of the impulse response of the channel in each polyphase branch
% the number of branches that have a coefficient less than the other
% branches is
n_short = mod(4-N, 4);
% then the number of coefficients in each branch is
N_i(1:4-n_short) = ceil(N/4);
N_i(4-n_short + 1:4) = ceil(N/4) - 1;


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
clear l r

x = [p; p(1:max(N_i)-1)];
x(x == 0) = -1;


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
Tc = 1;     % This is the smallest time interval we want to simulate
T = 4*Tc;   % Time sampling interval of the input of the channel
Kdb = 3;    % 3 dB
K = 10^(Kdb/10);
fd = 5*10^-3/T; % doppler frequency
Tq = Tc;    % Fundamental sampling time. This is the same as Tc, right???
% We stick to anastasopolous chugg paper (1997) and choose Tp such that
% fd*Tp = 0.1
tau_rms = 0.3*T;  % Tc = 1, as it is the fundamental sampling time
tau = 0:Tc:N_h-1;
Tp = 1/10 * (1/fd) * Tq; % Sampling time used for filtering the white noise
pdp_gauss = 1/tau_rms * exp(-tau/tau_rms);

C = sqrt(K/(K+1));
% normalize pdp: it must be sum(E[|gtilde_i|^2]) = 1 - C^2
pdp_gauss = pdp_gauss.*(1-C^2)/sum(pdp_gauss);

snr = 10; % db
snr_lin = 10^(snr/10);
sigma_w = 1/(T/Tc*snr); % the PN sequence has power 1


fprintf('fd * Tp = %d \n', Tp*fd);
g_samples_needed = 200000; % Some will be dropped because of transient
w_samples_needed = ceil(g_samples_needed / Tp);
transient = ceil(g_samples_needed/4);

% Generate complex-valued Gaussian white noise with zero mean and unit variance
%rng('default');
g_mat = zeros(N_h, g_samples_needed - transient);
for ray = 1:N_h
    w = wgn(w_samples_needed,1,0,'complex');
    %fprintf('variance of white noise=%d \n', var(w))
    
    % Filter the wgn with a narrowband filter. The filter will have the
    % classical Doppler spectrum in the frequency domain, with f_d * Tc = 1.25*10^-3
    % By using the approach suggested in anastchugg97 we use an iir filter
    % with known coefficients which is the convolution of a Cheby lowpass and a shaping filter.
    % Note that it is
    % possible to initialize properly the filters, consider doing it in
    % following versions since it allows to start in steady state conditions
    % and avoid dropping about many samples (it is a very cool thing!)
    gprime = filter(b_dopp, a_dopp, w);
    fprintf('g%d after interpolation has mean %d and variance %d  \n', ray, mean(gprime), var(gprime));
    %Gpr = fft(gprime);
    %figure, plot(20*log10(abs(Gpr)))
    
    % Interpolation
    t = 1:length(gprime);
    t_fine = Tq/Tp:Tq/Tp:length(gprime);
    
    g_fine = interp1(t, gprime, t_fine, 'spline');
    fprintf('g%d after interpolation has mean %d and variance %d  \n', ray, mean(g_fine), var(g_fine));
    % Drop the transient
    g_mat(ray, :) = g_fine(transient+1:end);
end

% Energy scaling
for k = 1:N_h
    g_mat(k, :) = g_mat(k, :)*sqrt(pdp_gauss(k));
end

% Only for LOS component
g_mat(1, :) = g_mat(1, :) + C;

clear tau tau_rms a_dopp b_dopp g_fine gprime h_dopp pdp_gauss t_dopp t_fine g_samples_needed w_samples_needed M_d k ray

%% Generate desired signal via a polyphase implementation
%!!! This implementation works for N_h <= 4 !!!

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
g_used_coeff = g_used_coeff(:, max(N_i) - 1 + 1:end); % N-1 is the max transient, +1 because of MATLAB

% compute the mean coefficient of impulse response of each ray over the interval
% of interest of L samples in order to compare them with the estimated
% impulse response
h_mean = mean(g_used_coeff, 2);

% create four different d_i vector, by sampling at step 4 the complete vector
% d. Each of them is the output of the polyphase brach at "lag" i
d_poly = zeros(length(d)/4, 4); % each column is a d_i
for idx  = 1:4
    d_poly(:, idx) = d(idx:4:end);
end

figure, stem(0:T:(length(d)-1), abs(x)), hold on, stem(0:Tc:length(d)-1, abs(d));
legend('x', 'd');

%% Estimate h

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

% tentative of computing d_hat
% h_hat_array = reshape(h_hat, 4*max(N_i), 1); 
% 
% d_hat = zeros(T*length(x) - N,1);
% h_used_coeff = zeros(4, length(x)); % check these dimensions
% for k = max(N_i):length(x)-1
%     % Generate white noise
%     for idx = 0:4-1 % actually idx varies between 0 and 4 since there are four branches
%         d_hat(k*T + idx*Tc + 1) = h_hat_array(idx*Tc+1:4:end).' * flipud(x(k+1-max(N_i) + 1:k+1));
%         % store the coefficient actually used, it will be useful later on
%         %h_used_coeff(idx + 1, k + 1) = g_mat(idx+1,k*T + idx*Tc+1);
%         
%     end
% end

%% Compute d_hat

% We need to discard N-1-(T/Tc-1) = N - T/Tc = N-4 samples for the
% transient. Actually we need to discard N_tr = max(0, N-4).
% %We already disregarded ?? floor((N-4)/4)
% x has L+max(N_i)-1 samples, we are considering L+1, that is we are
% disregarding max(N_i)-2 samples of x. This is equivalent to discarding
% 4*(max(N_i)-2) samples of d. We need to discard instead N-4 samples, so
% we still have to discard some other samples of d. How many?
% N_tr-4*max(N_i)+8 =
%   = N-4*(max(N_i)-1) if N>4
%   = -4*(max(N_i)-1) = 0  otherwise
% that practically is computed as max(0, N-4)-4*ceil(N/4)+8. Note that this
% yields a result that has a periodic ramp behaviour from 1 to 4, except
% for the first 4 values (N<=4) in which it is constantly 4.
% If max(N_i)=1 then we must keep all samples of x. What happens is that x
% has L samples, we take these L samples with a 0 in the front, h_hat is a
% column vector and x_toep is a row vector of length L+1. The resulting
% d_hat is as follows. The first column is all zeros and in this case is
% completely useless. The second column needs to be kept, since it is the
% output of the system in the first 4 time instants, and they are all
% useful because each branch of the filter has order 0 hence it does not
% depend on past values of x. Indeed, 4 is the value we get from the
% expression we derived above.
x = [p; p(1:max(N_i)-1)];
x(x == 0) = -1;
x = [0; x];
x_toep = toeplitz(x);
x_toep = x_toep(1:max(N_i), end-L:end);
% d_hat with the final part of the transient or with some useless zeros.
d_hat = h_hat * x_toep;
% Get d_hat in a line and then discard samples.
d_hat = reshape(d_hat, numel(d_hat), 1);
d_hat_discard_num = max(0, N-4)-4*ceil(N/4)+8;
d_hat = d_hat(d_hat_discard_num + 1 : end);
d_no_trans = d(end-length(d_hat)+1 : end);
% Plot
figure
subplot 211, plot([real(d_no_trans), real(d_hat)])
legend('Re[d]', 'Re[d_{hat}]'), grid on
title('Real part of actual and estimated desired signal')
subplot 212, plot([imag(d_no_trans), imag(d_hat)])
legend('Im[d]', 'Im[d_{hat}]'), grid on
title('Imaginary part of actual and estimated desired signal')


%% Error estimate

error_func = sum(abs(d_hat - d_no_trans).^2);
fprintf('EPSILON = %d\n', error_func);


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
