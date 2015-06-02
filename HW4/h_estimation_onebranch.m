function [ h_hat, d_hat ] = h_estimation_onebranch( x, d, L, N )
% This function performs the estimation of the h coefficients, given the
% input sequence, the output of the channel and the number of coefficients
% to estimate. Additionally, the function outputs d_hat, i.e. the output of
% a channel that would have the estimated coefficients as impulse response.

% In this case the estimation cannot be performed.
if N > L
    h_hat = [];
    d_hat = [];
    return
end

%% Estimate h


% Using the data matrix (page 246), easier implementation
h_hat = zeros(1,N);
I = zeros(L,N);
for column = 1:N
    I(:,column) = x(N-column+1:(N+L-column));
end
o = d(N:N + L - 1);

% Compute the Phi matrix and the theta vector
Phi = I'*I;
theta = I'*o;

h_hat(1, 1:N) = Phi \ theta;

%% Compute d_hat

d_hat = conv(x, h_hat);
d_hat = d_hat(N : N+L-1);

end