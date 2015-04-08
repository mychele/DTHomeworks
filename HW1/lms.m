% Implement LMS algorithm

close all;
clear all;
clc;

%% Load "continuous PSD" signal
load('split_signal.mat', 'z_continuous');
z = z_continuous - mean(z_continuous);
K = length(z); % signal length
autoc_z = autocorrelation(z, round(K/5));

%% AR
% the knee is apparently at N = 3, however LMS doesn't converge for N = 3
% in the required number of iterations
% compute the vector of coefficients a
N = 2;
[a, sigma_w] = arModel(N, autoc_z);
[H, omega] = freqz(1, [1; a], K, 'whole');

%%
upper_limit = 399; %MATLAB requires indices from 1 to 401
c = zeros(N, upper_limit + 1); % init c vector, no info -> set to 0
% each column of this matrix is c(k), a vector with coefficients from 1 to
% N (since we are implementing the predictor)!
e = zeros(1, upper_limit);

mu = 0.42/(autoc_z(1)*N); % actually mu must be > 0 and < 2/(N r_z(0))
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
hold on
plot(1:upper_limit, 10*log10(sigma_w)*ones(1, upper_limit))
title('Error function at each iteration');

% Find the value of coefficients at instant k = 350 and the average of e
% over k \in [350 - 10, 350 + 10]

ind = 350;
win_side_len = 10;
win_len = 2*win_side_len + 1;
c_350 = c(:, ind);
e_350_av = 10*log10(sum(abs(e(ind-win_side_len:ind+win_side_len)).^2)/win_len);
