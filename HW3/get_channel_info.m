function [ h_i, est_sigmaw ] = get_channel_info( r, N1, N2 )
%GET_CHANNEL_INFO

L = 15;
Nseq = 10;
trainingsymbols = ts_generation(L, Nseq);


% --- Estimate impulse response h @T and compute estimated noise power
N = N1+N2+1;
x_for_ls = trainingsymbols(end - (L+N-1) + 1 : end);
% r is in T
d_for_ls = r(end - (L+N-1) + 1 - N1 : end - N1 );  % is -N1 correct?
[h_i, r_hat] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);
d_no_trans = d_for_ls(N : N+L-1);
est_sigmaw = sum(abs(r_hat - d_no_trans).^2)/length(r_hat);

h_i = h_i.'; % for convenience

end