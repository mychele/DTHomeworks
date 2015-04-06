%% This is the implementation of the Recursive Least Squares algorithm (RLS)

% For reference, see pages 197, 201-203 of the Benvenuto-Cherubini book.

% Clear stuff
close all
clear all

%% Load data

% Load actual signal
z = load('data for hw1.mat');
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
z = z/30 + 3*exp(1i*2*pi*0.1*(1:1000).' + 0.1);
K = length(z); % signal length
autoc_z = autocorrelation(z, K/5);


%% Initialisation
N = 2;
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

%% Begin iterating
% Remember, we are implementing a predictor, so the z(k) of the book is 
% actually z(k-1) for us. See page 201 for reference.
%
% NOTE: I _hate_ MATLAB's indexing from 1. All indices are kept just like
% they are in the book, and k simply starts from 2 instead of 1. 

%c(:, 1) = [(15 + 2i); (2 -5i)];
w0 = 2*pi*0.1;
const = 2;
x = (const * exp(1i * w0 * (1 : upper_limit+1))).';

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
    y = x(k) * (c(1, k-1) - 1i * c(2, k-1));
    % Compute a priori estimation error (with old coefficients)
    epsilon(k) = d(k) - y;
    
    %epsilon(k) = d(k) - z_k.' * c(:, k-1);
    c(:, k) = c(:, k-1) + epsilon(k) * k_star(:,k);
    %e(k) = d(k) - x_k.' * c(:,k);
    
    % Output y(k) computed with new coefficients c(k)
    y = x(k) * (c(1, k) - 1i * c(2, k));
    % Compute a posteriori estimation error (with new coefficients)
    e(k) = d(k) - y;
    
    P(:,:,k) = 1/lambda * (P(:,:,k-1) - k_star(:,k)*pi_star(:,k)');
end

% End of computation.

% Plot c coefficients
for index = 1:N
    figure
    subplot(2, 1, 1)
    plot(1:upper_limit+1, real(c(index, :)))
    title(['Real part of c' int2str(index)]);
    subplot(2, 1, 2)
    plot(1:upper_limit+1, imag(c(index, :)))
    title(['Imaginary part of c' int2str(index)]);
end

% Plot the error.
figure, plot(1:upper_limit+1, 10*log10(abs(e).^2))
title('Error function at each iteration')
xlabel('Iteration (k)'), ylabel('|e(k)|^2 (db)')

figure
plot(1:upper_limit+1, 20*log10((abs(e))' ./ abs(d)))
title('Ratio between |e(k)|^2 and |z(k)|^2')
xlabel('Iteration (k)'), ylabel('Ratio (dB)')


%% Find amp and phase

% Average of coefficients from some iteration on, when hopefully they have converged
expcoeff = mean(c(:, floor(upper_limit*0.9) : upper_limit), 2)

% Acos(t) + jA*sin(t)
% = c1 * cos(t0) + c2 * sin(t0) + j (c1 * cos(t0 - pi/2) + c2 * sin(t0 - pi/2))
% = c1 * (cos(t0) + j * sin(t0)) + c2 * (sin(t0) - j * cos(t0))
% = exp(j * t0) (c1 - j*c2)
% Thus, Acos(t) = c1 cos(t0) + c2 sin(t0)
estimatedsine = x * (expcoeff(1) - 1i*expcoeff(2));
amp = const * abs(expcoeff(1) - 1i*expcoeff(2))


figure, plot(real(estimatedsine)), hold on, plot(imag(estimatedsine), 'r')
title('Estimated signal - imag and real parts')
legend('Real part', 'Imag part')
figure, plot3(1:30, real(estimatedsine(1:30)), imag(estimatedsine(1:30)))
title('First 30 samples of estimated signal')

figure
subplot(2, 1, 1), plot(real(estimatedsine)), hold on, plot(real(z), 'r')
title('Real part of estimated exp vs total signal')
legend('Estimated exp', 'Total signal')
subplot(2, 1, 2), plot(imag(estimatedsine)), hold on, plot(imag(z), 'r')
title('Imag part of estimated exp vs total signal')
legend('Estimated exp', 'Total signal')