% Implement LMS algorithm

close all
clear all
clc
rng default

ctot = zeros(3, 1000);

iterations = 200;
for i=1:iterations
    
    %% Load data
    z = randn(5000, 1);
    z = filter(1, [1 0.2-0.2i 0.5 0.2], z);
    K = length(z); % signal length
    autoc_z = autocorrelation(z, K/10);
    
    %% AR
    % the knee is apparently at N = 3
    % compute the vector of coefficients a
    N = 3;
    [a, sigma_w] = arModel(N, autoc_z);
    %[H, omega] = freqz(1, [1; a], K, 'whole');
    %hold on, plot(omega*10000/2/pi, 10*log10(abs(H)))
    
    %%
    N = 3; % order of the predictor
    upper_limit = 999; %MATLAB requires indices from 1 to 401
    c = zeros(N, upper_limit + 1); % init c vector, no info -> set to 0
    % each column of this matrix is c(k), a vector with coefficients from 1 to N (since we are implementing the predictor)!
    e = zeros(1, upper_limit);
    
    mu = 0.1/(autoc_z(1)*N); % actually mu must be > 0 and < 2/(N r_z(0))
    
    % watch out, in the predictor y(k) = transp(x(k-1))c(k)
    for k = 1:upper_limit
        if (k < N + 1)
            z_k_1 = flipud([zeros(N - k + 1, 1); z(1:k - 1)]); % input vector z_vec_(k-1) of length N
            % for k = 1 z(1:0) is an empty matrix
        else
            z_k_1 = flipud(z((k - N):(k-1))); % we need the input from k - 1 to k - N
        end
        y_k = z_k_1.'*c(:, k);
        e_k = z(k) - y_k; % the reference signal d(k) is actually the input at sample k
        e(k) = e_k;
        c(:, k + 1) = c(:, k) + mu*e_k*conj(z_k_1); % update the filter, c(k+1) = c(k) + mu*e(k)*conj(z(k-1))
    end
    ctot = ctot + c / iterations;
    
    disp(i)
end



figure
subplot(2, 1, 1)
plot(1:upper_limit+1, real(ctot(1, :)), [1, upper_limit+1], -real(a(1))* [1 1])
title('Real part of c1');
subplot(2, 1, 2)
plot(1:upper_limit+1, imag(ctot(1, :)), 1:upper_limit+1, -imag(a(1)))
title('Imaginary part of c1');

figure
subplot(2, 1, 1)
plot(1:upper_limit+1, real(ctot(2, :)), 1:upper_limit+1, -real(a(2)))
title('Real part of c2');
subplot(2, 1, 2)
plot(1:upper_limit+1, imag(ctot(2, :)), 1:upper_limit+1, -imag(a(2)))
title('Imaginary part of c2');
figure
subplot(2, 1, 1)
plot(1:upper_limit+1, real(ctot(3, :)), 1:upper_limit+1, -real(a(3)))
title('Real part of c3');
subplot(2, 1, 2)
plot(1:upper_limit+1, imag(ctot(3, :)), 1:upper_limit+1, -imag(a(3)))
title('Imaginary part of c3');


return


figure
subplot(2, 1, 1)
plot(1:upper_limit+1, real(c(1, :)), [1, upper_limit+1], -real(a(1))* [1 1])
title('Real part of c1');
subplot(2, 1, 2)
plot(1:upper_limit+1, imag(c(1, :)), 1:upper_limit+1, -imag(a(1)))
title('Imaginary part of c1');

figure
subplot(2, 1, 1)
plot(1:upper_limit+1, real(c(2, :)), 1:upper_limit+1, -real(a(2)))
title('Real part of c2');
subplot(2, 1, 2)
plot(1:upper_limit+1, imag(c(2, :)), 1:upper_limit+1, -imag(a(2)))
title('Imaginary part of c2');
figure
subplot(2, 1, 1)
plot(1:upper_limit+1, real(c(3, :)), 1:upper_limit+1, -real(a(3)))
title('Real part of c3');
subplot(2, 1, 2)
plot(1:upper_limit+1, imag(c(3, :)), 1:upper_limit+1, -imag(a(3)))
title('Imaginary part of c3');

figure, plot(1:upper_limit, 10*log10(abs(e).^2))
title('Error function at each iteration');

