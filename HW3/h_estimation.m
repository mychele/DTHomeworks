function [ h_hat, d_hat ] = h_estimation( x, d, L, N_i )
% This function performs the estimation of the h coefficients, given the
% input sequence, the output of the channel and the number of coefficients 
% to estimate. Additionally, the function outputs d_hat, i.e. the output of
% a channel that would have the estimated coefficients as impulse response.

% In this case the estimation cannot be performed.
if max(N_i) > L
    h_hat = [];
    d_hat = [];
    return
end

%% Estimate h

% Create four different d_i vectors, by sampling with step 4 the complete
% vector d. Each of them is the output of the branch that has "lag" iTc.
d_poly = zeros(ceil(length(d)/4), 4); % each column is a d_i
for idx  = 1:4
    d_poly(:, idx) = d(idx:4:end);
end

% Using the data matrix (page 246), easier implementation
h_hat = zeros(4,max(N_i));
% We're estimating 4 branches of the polyphase representation.
% This matrix has the maximum number of coefficients for each of the four
% branches. The unused (i.e. unestimated) ones will be left zero.
for idx = 1:4
    if N_i(idx) > 0
        I = zeros(L,N_i(idx));
        for column = 1:N_i(idx)
            I(:,column) = x(N_i(idx)-column+1:(N_i(idx)+L-column));
        end
        o = d_poly(N_i(idx):N_i(idx) + L - 1, idx);
        
        % Compute the Phi matrix and the theta vector
        Phi = I'*I;
        theta = I'*o;
        
        h_hat(idx, 1:N_i(idx)) = Phi \ theta;
    end % if N_branch is 0 don't estimate and leave hhat to 0
end

%% Compute d_hat

x = [0; x];
x_toep = toeplitz(x);
x_toep = x_toep(1:max(N_i), end-L:end);
% d_hat with the final part of the transient or with some useless zeros.
d_hat = h_hat * x_toep;
% Get d_hat in a line and then discard samples.
d_hat = reshape(d_hat, numel(d_hat), 1);
d_hat_discard_num = max(0, sum(N_i)-4)-4*ceil(sum(N_i)/4)+8;
d_hat = d_hat(d_hat_discard_num + 1 : end);
d_no_trans = d(end-length(d_hat)+1 : end);

end