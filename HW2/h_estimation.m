function [ h_hat, d_hat ] = h_estimation( x, d, L, N_i )
%H_ESTIMATION

% In this case the estimation cannot be performed.
if max(N_i) > L
    h_hat = [];
    d_hat = [];
    return
end

%% Estimate h

% Create four different d_i vectors, by sampling with step 4 the complete
% vector d. Each of them is the output of the branch that has "lag" i.
d_poly = zeros(length(d)/4, 4); % each column is a d_i
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
        
        Phi = I'*I;
        theta = I'*o;
        
        h_hat(idx, 1:N_i(idx)) = Phi \ theta;
    end % if N_branch is 0 don't estimate and leave hhat to 0
end


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




% Plot (debug purposes)
% figure
% subplot 211, plot([real(d_no_trans), real(d_hat)])
% legend('Re[d]', 'Re[d_{hat}]'), grid on
% title('Real part of actual and estimated desired signal')
% subplot 212, plot([imag(d_no_trans), imag(d_hat)])
% legend('Im[d]', 'Im[d_{hat}]'), grid on
% title('Imaginary part of actual and estimated desired signal')


end