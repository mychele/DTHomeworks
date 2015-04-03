% Implement LMS algorithm

close all;
clear all;
clc;

%% Load data
z = load('data for hw1.mat');
fir_bs_1 = load('fir_bs_1.mat');
firbs = fir_bs_1.fir_bs_1;
z = z.z.'; % make a column vector
z = z - mean(z); % remove average
K = length(z); % signal length
autoc_z = autocorrelation(z, K/5);

%% AR
% the knee is apparently at N = 3
% compute the vector of coefficients a
N = 2;
[a, sigma_w] = arModel(N, autoc_z);
[H, omega] = freqz(1, [1; a], K, 'whole');

%%
upper_limit = 999; %MATLAB requires indices from 1 to 401
c = zeros(N, upper_limit + 1); % init c vector, no info -> set to 0
%c(:,1) = -a +5;
% each column of this matrix is c(k), a vector with coefficients from 1 to
% N (since we are implementing the predictor)!
e = zeros(1, upper_limit);


mu = 0.1/(autoc_z(1)*N); % actually mu must be > 0 and < 2/(N r_z(0))

% watch out, in the predictor y(k) = transp(x(k-1))c(k)
for k = 1:upper_limit
    if (k < N + 1)
        z_k_1 = flipud([zeros(N - k + 1, 1); z(1:k - 1)]); % input vector z_vec_(k-1) of length N
        % for k = 1 z(1:0) is an empty matrix
        y_k = z_k_1.'*c(:, k);
    else
        z_k_1 = flipud(z((k - N):(k-1))); % we need the input from k - 1 to k - N
        y_k = z_k_1.'*c(:, k);
    end
    e_k = z(k) - y_k; % the reference signal d(k) is actually the input at sample k
    e(k) = e_k;
    c(:, k + 1) = c(:, k) + mu*e_k*conj(z_k_1); % update the filter, c(k+1) = c(k) + mu*e(k)*conj(z(k-1))
end

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


figure, plot(1:upper_limit, 10*log10(abs(e).^2))
title('Error function at each iteration');

