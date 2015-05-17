function [ h_i, est_sigmaw, N1, N2] = get_channel_info( r, N1, N2, assumed_m_opt, T )
%GET_CHANNEL_INFO

L = 15;
Nseq = 10;
trainingsymbols = ts_generation(L, Nseq);


% % --- Estimate optimal timing phase
% m_min = 0;
% m_max = 43;
% i = 1;
% crossvec = zeros(m_max - m_min, 1);
% for m = m_min:m_max
%     r_part = r(m+1:T:m+1+T*(L-1)); % pick L samples from r, spaced by T
%     crossvec(i) = abs(sum(r_part.*conj(trainingsymbols(1:L)))/L); % as in 7.269
%     i = i + 1;
% end
% [~, m_opt] = max(crossvec);
% m_opt = m_opt - 1; % because of MATLAB indexing
% init_offs = mod(m_opt, T);
% h_center = floor(m_opt/ T);
% N1 = h_center;
% N2 = length(h_i) - N1 - 1;


% --- Estimate impulse response h @T and compute estimated noise power

t0 = floor(assumed_m_opt / T);
init_offs = mod(assumed_m_opt, T);
N = N1+N2+1;
x_for_ls = trainingsymbols(end-(L+N-1) + 1 : end);
d_for_ls = r(end - T*(L+N-1) + 1 + init_offs - T*N1: T :end - T*N1);
[h_i, r_hat] = h_estimation_onebranch(x_for_ls, d_for_ls, L, N);
d_no_trans = d_for_ls(N : N+L-1);
est_sigmaw = sum(abs(r_hat - d_no_trans).^2)/length(r_hat);

h_i = h_i.'; % for convenience

end