%function [decisions] = DFE(x, h, N1, N2)
% Arguments:
%   x: input sequence
%   h: channel impulse response (actually, its estimate, h_hat)
%   N1: number of precursors
%   N2: number of postcursors

receiver_util;

snr_col = 1;

rT = r_T4(init_offs+1:4:end, snr_col); % data sampled in T
x = rT/h(N1+1).'; % data normalized by h0
hi = h/h(N1+1).'; % impulse response normalized by h0

% Initialize useful quantities
sigma_a = 2;
N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1 = N1+N2+1;   % FF filter: equal to the span of h
M2 = M1-1;      % FB filter: one less than the FF filter
D = (N-1)/2; %(M1-1);   % D is chosen large first and then decreased
K = length(x);
a_k = zeros(K,1);

% Zero padding of the i.r.
nb0 = 20;
nf0 = 20;
hi = [zeros(nb0,1); hi; zeros(nf0,1)];

% Get the Weiner-Hopf solution
p = zeros(M1, 1);
for i = 0:(M1 - 1)
    p(i+1) = sigma_a * conj(hi(N1+nb0+1+D-i));
end

R = zeros(M1);
for row = 0:(M1-1)
    for col = 0:(M1-1)
        first_sum = (hi((nb0+1):(N1+N2+nb0+1))).' * ...
            conj(hi((nb0+1-(row-col)):(N1+N2+nb0+1-(row-col))));
        second_sum = (hi((N1+nb0+1+1+D-col):(N1+nb0+1+M2+D-col))).' * ...
            conj((hi((N1+nb0+1+1+D-row):(N1+nb0+1+M2+D-row))));
        r_w = (row == col) * est_sigmaw; % This is a delta only if there is no g_M.
        
        R(row+1, col+1) = sigma_a * (first_sum - second_sum) + r_w;
        
    end
end

c_opt = R \ p;

b = zeros(M2,1);
for i = 1:M2
    b(i) = - (fliplr(c_opt.')*hi((i+D+N1+nb0+1-M1+1):(i+D+N1+nb0+1)));
end

%% TODO plot hhat, c, psi=conv(h,c), b and get a sense of what is happening
figure
subplot(4,1,1)
stem(abs(h))
title('|h_{hat}|')
xlim([1 14])
ylim([0 1])
subplot(4,1,2)
stem(abs(c_opt))
title('|c|')
xlim([1 14])
ylim([0 1])
subplot(4,1,3)
stem(abs(conv(h,c_opt)))
title('|psi|')
xlim([1 14])
ylim([0 1])
subplot(4,1,4)
stem(abs(b))
title('|b|')
xlim([1 14])
ylim([0 1])

%% Threshold detector

% known data
preamble = MLsequence(15);

y = zeros(length(x),1);
detected = zeros(length(x), 1);
for k = 0:length(x)-1
    if (k < M1 - 1)
        xconv = [flipud(x(1:k+1)); zeros(M1 - k - 1, 1)];
    else
        xconv = flipud(x(k-M1+1 + 1:k + 1));
    end
    
    if (k <= M2)
        a_old = [flipud(detected(1:k)); zeros(M2 - k, 1)];
    else
        a_old = flipud(detected(k-M2+1:k));
    end
    
    y(k+1) = c_opt.'*xconv + b.'*a_old;
    detected(k+1) = qpsk_td(y(k+1));
end

figure
subplot(4,1,1)
stem(packet)
xlim([0,length(packet)])
title('Training Symbols')
subplot(4,1,2)
stem(x(t0+1:end))
xlim([0,length(packet)])
title('Received Symbols')
subplot(4,1,3)
stem(y(t0+D+1:end))
xlim([0,length(packet)])
title('Filtered Symbols')
subplot(4,1,4)
stem(detected(t0+D+1:end))
xlim([0,length(packet)])
title('Detected Symbols')

% figure
% for i = D + t0 + 1:length(trainingsymbols)
%     plot(trainingsymbols(i-D-t0), 'ok'), hold on, plot(x(i-D), 'or'), 
%     hold on, plot(y(i), 'og'),%  plot(detected(i), 'ob')
%     xlim([-2, 2]), ylim([-2, 2]), grid on;
%     legend('sent', 'received', 'filtered', 'detected')
%     pause
%     hold off
% end

% TODO check that the R matrix should be Hermitian and Toeplitz for a
% LE, while it should be only Hermitian for a DFE
% Threshold to take the decision
decisions = detected(t0+D+1:end);
%end