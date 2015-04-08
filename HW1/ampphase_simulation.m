%% This is a simulation designed to test how well does the RLS method works 
% when it comes to find amplitude, phase of a spectral line, given an
% approximated frequency derived from the observation of the peak in the
% DFT. The idea is to apply the computation for frequencies in an interval
% around the approximated one and pick the one that gives the best
% crosscorrleation(0) between the signal we created and the reconstructed
% signal. We will see that for the freq which gives the best xcorr(0) the
% estimates of given amplitude and phase are correct. 

% For reference, see pages 197, 201-203 of the Benvenuto-Cherubini book.

% Show that 1 coefficient is enough
% Ae^(j w0 k + j phi) = Ae^jphi e^jw0k, c = Ae^jphi

% Alternative, long and useless proof
% Ae^j(w0k + phi)
% Acos(w0k + phi) + jAsin(w0k + phi)

% Acos(w0k)cos(phi) - Asin(w0K)sin(phi) + jAsin(w0k)cos(phi) +
% jAcos(w0k)sin(phi)

% cm = Acos(phi), cd = Asin(phi)
% cm^2 + cd^2 = A^2(cos^2 + sin^2) = A^2

%cmcos(w0k) + jcmsin(w0k) + jcdcos(w0k) - cdsin(w0k))

%cm e^(jw0k) + j cd e^(jw0k)

% (cm + j cd) e^(jw0k)


% Clear stuff

close all
clear all

sim_length = 1000;

% Set our parameters
amp_est = zeros(sim_length, 1);
phi_est = zeros(sim_length, 1);
rng('default'); % creates reproducible results
w0 = rand();
r_phi = 2*pi*rand();
r_amp = 10*rand();

% Simulate sim_length times
for index = 1:sim_length
    z = wgn(1000, 1, 10) + r_amp*exp(1i*2*pi*w0*(1:1000).' + 1i * r_phi);
    K = length(z);
    autoc_z = autocorrelation(z, K/5);
    corr_vec = zeros(6, 1);
    amp_vec = zeros(6, 1);
    phi_vec = zeros(6,1);
    i = 1;
    for w1 = w0-0.02:0.01:w0+0.02
        
        
        % Initialisation
        N = 1; % see the first comment
        upper_limit = length(z)-1;%399; % Number of iterations of the algorithm
        lambda = 1; % Forgetting factor. For 1, we do not forget past values
        c = zeros(N, upper_limit+1); % Coefficient vector
        delta = autoc_z(1)/100; % Value at which to initialise P
        % P is a N+1 square matrix. P(n) is achieved by making P a parallelogram
        % Access P by using P(row, column, time)
        P(:,:,1) = (1/delta) * eye(N);
        pi_star = zeros(N, upper_limit+1); % pi_star is a series of column vectors
        r = zeros(1,upper_limit+1);   % r is a vector of scalars
        k_star = zeros(N, upper_limit+1);
        d = z; % The reference signal is the input at time k
        epsilon = zeros(1, upper_limit+1); % The a priori estimation error
        e = zeros(1,upper_limit+1);
        
        % Begin iterating
        % Remember, we are implementing a predictor, so the z(k) of the book is
        % actually z(k-1) for us. See page 201 for reference.
        %
        % NOTE: I _hate_ MATLAB's indexing from 1. All indices are kept just like
        % they are in the book, and k simply starts from 2 instead of 1.
        
        w = 2*pi*w1;
        const = 1;
        x = (const * exp(1i * w * (1 : upper_limit+1))).';
        
        for k = 2:upper_limit+1
            % Cut off the x(k-1) for this iteration (this part is stolen from the
            % lms implementation), handling the case in which k < N.
            if (k < N)  % Fill up with zeros
                x_k = flipud([zeros(N - k, 1); x(1:k)]);
            else % Just cut the input vector
                x_k = flipud(x((k - N + 1):(k)));
            end
            pi_star(:,k) = P(:,:,k-1) * conj(x_k);
            r(k) = 1/(lambda + x_k.' * pi_star(:,k));
            k_star(:,k) = r(k) * pi_star(:,k);
            
            % Output y(k) computed with old coefficients c(k-1)
            y = x(k) * (c(1, k-1));
            % Compute a priori estimation error (with old coefficients)
            epsilon(k) = d(k) - y;
            
            c(:, k) = c(:, k-1) + epsilon(k) * k_star(:,k);
            
            % Output y(k) computed with new coefficients c(k)
            y = x(k) * (c(1, k));
            % Compute a posteriori estimation error (with new coefficients)
            e(k) = d(k) - y;
            
            P(:,:,k) = 1/lambda * (P(:,:,k-1) - k_star(:,k)*pi_star(:,k)');
        end
        
        % End of computation.
              
        % Find amp and phase
        
        % Average of coefficients from some iteration on, when hopefully they have converged
        expcoeff = mean(c(:, floor(upper_limit*0.9) : upper_limit), 2);
        
        estimatedsine = x * (expcoeff(1));
        corr = crosscorrelation(estimatedsine, z, length(z)/5);
        corr_vec(i) = corr(1);
        amp_vec(i) = const*abs(expcoeff(1));
        phi_vec(i) = angle(expcoeff(1));
        i = i+1;
    end
    [mx, j] = max(abs(corr_vec));
    amp_est(index) = amp_vec(j);
    phi_est(index) = phi_vec(j);
end

%% Statistical things

mse_amp = sum((amp_est-r_amp).^2)/length(amp_est);
% watch out for the following, if the r_phi is over pi then use (-2*pi+r_phi)
mse_phi = sum((phi_est-(-2*pi+r_phi)).^2)/length(amp_est);

%% Different approach
% Check if it picks the correct frequency. Freq, amp and phase will be different
% each time.
% The correct index that should appear in ind_j is span/step + 1 (the
% center of the vector of frequencies passed to RLS, which is actually the
% freq of the input signal)

sim_length = 1000;

% Set our parameters
amp_est_2 = zeros(sim_length, 1);
phi_est_2 = zeros(sim_length, 1);
ind_j = zeros(sim_length, 1);
rng('default');
span = 0.05;
step = 0.01;

% Simulate sim_length times
for index = 1:sim_length
    w0 = rand();
    r_phi = pi*rand() - pi; %[-pi, pi] phase
    r_amp = 10*rand();
    z = wgn(1000, 1, 10) + r_amp*exp(1i*2*pi*w0*(1:1000).' + 1i * r_phi);
    K = length(z);
    autoc_z = autocorrelation(z, K/5);
    
    corr_vec = zeros(2*(span/step) + 1, 1);
    amp_vec = zeros(2*(span/step) + 1, 1);
    phi_vec = zeros(2*(span/step) + 1, 1);
    i = 1;
    for w1 = w0-span:step:w0+span % Then w0 should be the 6th element of the vector (span/step + 1)
        
        % Initialisation
        N = 1; % see the first comment
        upper_limit = length(z)-1;%399; % Number of iterations of the algorithm
        lambda = 1; % Forgetting factor. For 1, we do not forget past values
        c = zeros(N, upper_limit+1); % Coefficient vector
        delta = autoc_z(1)/100; % Value at which to initialise P
        % P is a N+1 square matrix. P(n) is achieved by making P a parallelogram
        % Access P by using P(row, column, time)
        P(:,:,1) = (1/delta) * eye(N);
        pi_star = zeros(N, upper_limit+1); % pi_star is a series of column vectors
        r = zeros(1,upper_limit+1);   % r is a vector of scalars
        k_star = zeros(N, upper_limit+1);
        d = z; % The reference signal is the input at time k
        epsilon = zeros(1, upper_limit+1); % The a priori estimation error
        e = zeros(1,upper_limit+1);
        
        % Begin iterating
        % Remember, we are implementing a predictor, so the z(k) of the book is
        % actually z(k-1) for us. See page 201 for reference.
        %
        % NOTE: I _hate_ MATLAB's indexing from 1. All indices are kept just like
        % they are in the book, and k simply starts from 2 instead of 1.
        
        %c(:, 1) = 15 + 2i;
        w = 2*pi*w1;
        const = 1;
        x = (const * exp(1i * w * (1 : upper_limit+1))).';
        
        for k = 2:upper_limit+1
            % Cut off the x(k-1) for this iteration (this part is stolen from the
            % lms implementation), handling the case in which k < N.
            if (k < N)  % Fill up with zeros
                x_k = flipud([zeros(N - k, 1); x(1:k)]);
            else % Just cut the input vector
                x_k = flipud(x((k - N + 1):(k)));
            end
            pi_star(:,k) = P(:,:,k-1) * conj(x_k);
            r(k) = 1/(lambda + x_k.' * pi_star(:,k));
            k_star(:,k) = r(k) * pi_star(:,k);
            
            % Output y(k) computed with old coefficients c(k-1)
            y = x(k) * (c(1, k-1));
            % Compute a priori estimation error (with old coefficients)
            epsilon(k) = d(k) - y;
            
            c(:, k) = c(:, k-1) + epsilon(k) * k_star(:,k);
            
            % Output y(k) computed with new coefficients c(k)
            y = x(k) * (c(1, k));
            % Compute a posteriori estimation error (with new coefficients)
            e(k) = d(k) - y;
            
            P(:,:,k) = 1/lambda * (P(:,:,k-1) - k_star(:,k)*pi_star(:,k)');
        end
        
        % End of computation.
              
        % Find amp and phase
        
        % Average of coefficients from some iteration on, when hopefully they have converged
        expcoeff = mean(c(:, floor(upper_limit*0.9) : upper_limit), 2);
        
        estimatedsine = x * (expcoeff(1));
        corr = crosscorrelation(estimatedsine, z, length(z)/5);
        corr_vec(i) = corr(1);
        amp_vec(i) = const*abs(expcoeff(1));
        phi_vec(i) = angle(expcoeff(1));
        i = i+1;
    end
    [mx, j] = max(abs(corr_vec));
    ind_j(index) = j;
    amp_est_2(index) = amp_vec(j) - r_amp;
    phi_est_2(index) = phi_vec(j) - r_phi;
end

find(ind_j ~= span/step + 1);
