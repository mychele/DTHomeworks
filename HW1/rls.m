%% This is the implementation of the Recursive Least Squares algorithm (RLS)

% We are trying to estimate the vector of coefficients c by an LS method.
% As a ballpark figure, this should converge ~10 times faster than the LMS
% algorithm. This comes at the cost of increased computational complexity.

% For reference, see pages 197, 201-203 of the Benvenuto-Cherubini book.

% Clear stuff
close all;
clear all;
clc;

%% Load "continuous PSD" signal
load('split_signal.mat', 'z_continuous');
z = z_continuous - mean(z_continuous);
K = length(z); % signal length
autoc_z = autocorrelation(z, round(K/5));

% Uncomment to load Dittadi's filtered white noise 
% z = randn(5000, 1);
% filtercoeff = [1, 0.2-0.5i, 0.2, 0.2];
% z = filter(1, filtercoeff, z);
% K = length(z); % signal length
% autoc_z = autocorrelation(z, K/10);


%% AR model
% compute the vector of coefficients a
N = 2;
[a, sigma_w] = arModel(N, autoc_z);
[H, omega] = freqz(1, [1; a], K, 'whole');


%% Initialisation
upper_limit = 399; % Number of iterations of the algorithm
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
epsilon = zeros(1, upper_limit+1); % The a posteriori estimation error
e = zeros(1,upper_limit+1);

%% Begin iterating
% Remember, we are implementing a predictor, so the z(k) of the book is 
% actually z(k-1) for us. See page 201 for reference.
%
% NOTE: I _hate_ MATLAB's indexing from 1. All indices are kept just like
% they are in the book, and k simply starts from 2 instead of 1. 
for k = 2:upper_limit+1
    % Cut off the x(k-1) for this iteration (this part is stolen from the 
    % lms implementation), handling the case in which k < N + 1.
    if (k < N + 1)  % Fill up with zeros
        z_k_1 = flipud([zeros(N - k + 1, 1); z(1:k - 1)]);
    else % Just cut the input vector
        z_k_1 = flipud(z((k - N):(k-1)));
    end
    pi_star(:,k) = P(:,:,k-1) * conj(z_k_1); 
    r(k) = 1/(lambda + z_k_1.' * pi_star(:,k));
    k_star(:,k) = r(k) * pi_star(:,k);
    epsilon(k) = d(k) - z_k_1.' * c(:, k-1);
    c(:, k) = c(:, k-1) + epsilon(k) * k_star(:,k);
    e(k) = d(k) - z_k_1.' * c(:,k);
    P(:,:,k) = 1/lambda * (P(:,:,k-1) - k_star(:,k)*pi_star(:,k)');
end

% End of computation.

% Plot c coefficients
for index = 1:N
    figure
    subplot(2, 1, 1)
    plot(1:upper_limit+1, real(c(index, :)), [1, upper_limit+1], -real(a(index))* [1 1])
    title(['Real part of c' int2str(index)]);
    subplot(2, 1, 2)
    plot(1:upper_limit+1, imag(c(index, :)), [1, upper_limit+1], -imag(a(index))* [1 1])
    title(['Imaginary part of c' int2str(index)]);
end

% Plot the error.
figure, plot(1:upper_limit+1, 10*log10(abs(e).^2), [1, upper_limit+1], 10*log10(sigma_w)*[1 1])
title('Error function at each iteration');

