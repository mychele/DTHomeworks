%% Amp-Phase Estimation RLS
% For reference, see pages 197, 201-203 of the Benvenuto-Cherubini book.

% Clear stuff
close all
clear all

%% Load data

% Load spectral signal
% load('split_signal.mat');
% z = z_lines;
% K = length(z); % signal length
% autoc_z = autocorrelation(z, floor(K/5));

% Load complete signal
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
K = length(z); % signal length
autoc_z = autocorrelation(z, floor(K/5));

f0 = 0.77; %estimated freq from DFT
w0 = 2*pi*f0;

span = 0.005;
step = 0.0001;

% vectors that will host temporary results
corr_vec = zeros(2*(span/step) + 1, 1);
amp_vec = zeros(2*(span/step) + 1, 1);
phi_vec = zeros(2*(span/step) + 1, 1);
i = 1;

for f1 = (f0 - span):step:(f0 + span)
    
    
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
    % NOTE: All indices are kept just like they are in the book, and k 
    % simply starts from 2 instead of 1.
    
    w = 2*pi*f1;
    const = 1;
    x = (const * exp(1i * w * (1 : upper_limit+1))).'; % reference signal
    
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
    corr = crosscorrelation(estimatedsine, z, floor(length(z)/5));
    corr_vec(i) = corr(1);
    amp_vec(i) = const*abs(expcoeff(1));
    phi_vec(i) = angle(expcoeff(1));
    i = i+1;
end
[mx, j] = max(abs(corr_vec));
amp_est = amp_vec(j);
phi_est = phi_vec(j);

% Plotting
figure, plot(real(estimatedsine)), hold on, plot(imag(estimatedsine), 'r')
title('Estimated signal - imag and real parts')
legend('Real part', 'Imag part')
figure, plot3(1:30, real(estimatedsine(1:30)), imag(estimatedsine(1:30)))
title('First 30 samples of estimated signal')