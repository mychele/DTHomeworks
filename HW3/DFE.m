%function [decisions] = DFE(x, h, N1, N2)
% Arguments:
%   x: input sequence
%   h: channel impulse response (actually, its estimate, h_hat)
%   N1: number of precursors
%   N2: number of postcursors

clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr = 10; % 6, 8, 10, 12, 14 % dB
L_data = 127;

%% Create, send and receive data with the given channel
[packet, r_T4, ~] = txrc(L_data, snr, T, Tc);

% estimate the channel using the first 100 samples (4*length(ts))
N = 7;
[ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);

% sample to get r @ T
init_offs = mod(m_opt, 4); % in T/4
t0 = floor(m_opt/4); % from now consider T = 1
T = 1;

%% Detection begins

rT = r_T4(init_offs+1:4:end); % data sampled in T
x = rT/h(N1+1).'; % data normalized by h0
hi = h/h(N1+1).'; % impulse response normalized by h0

% Initialize useful quantities
sigma_a = 2;
N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1 = N1+N2+1;   % FF filter: equal to the span of h
M2 = 0;      % FB filter: one less than the FF filter
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
len_plot = 14;
figure
subplot(4,1,1)
stem(-N1:N2, abs(hi(nb0+1:nb0+1+N1+N2)))
title('|h_{hat} normalized|')
xlim([-5 8])
ylim([0 1])
subplot(4,1,2)
stem(abs(c_opt))
title('|c|')
xlim([1 14])
ylim([0 1])
subplot(4,1,3)
stem(abs(conv(hi,c_opt)))
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

decisions = detected(t0+D+1:end-1);


figure
subplot(4,1,1)
stem(real(packet))
xlim([0,length(packet)])
title('Sent Symbols (real part)')
subplot(4,1,2)
stem(real(x(t0+1:end)))
xlim([0,length(packet)])
title('Received Symbols (real part)')
subplot(4,1,3)
stem(real(y(t0+D+1:end)))
xlim([0,length(packet)])
title('Filtered Symbols (real part)')
subplot(4,1,4)
stem(real(detected(t0+D+1:end))), hold on
stem(real(packet))
xlim([0,length(decisions)])
legend('Detected', 'Real')
title('Detected Symbols (real part)')

figure
subplot(4,1,1)
stem(imag(packet))
xlim([0,length(packet)])
title('Sent Symbols (imag part)')
subplot(4,1,2)
stem(imag(x(t0+1:end)))
xlim([0,length(packet)])
title('Received Symbols (imag part)')
subplot(4,1,3)
stem(imag(y(t0+D+1:end)))
xlim([0,length(packet)])
title('Filtered Symbols (imag part)')
subplot(4,1,4)
stem(imag(detected(t0+D+1:end))), hold on,
stem(imag(packet))
xlim([0,length(decisions)])
legend('Detected', 'Real')
title('Detected Symbols (imag part)')

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
%end

sent_bit = ibmap(packet);
rec_bit = ibmap(decisions);

figure, 
stem(sent_bit), hold on, stem(rec_bit)
legend('sent bit', 'received bit')



pbit = BER(packet, decisions);

