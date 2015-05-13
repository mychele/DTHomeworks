%function [decisions] = DFE(x, h, N1, N2)
% Arguments:
%   x: input sequence
%   h: channel impulse response (actually, its estimate, h_hat)
%   N1: number of precursors
%   N2: number of postcursors

rT = r(init_offs+1:4:end);
x = rT;
h = hi;

% Initialize useful quantities
sigma_a = 2;
N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
M1 = N1+N2+1;   % FF filter: equal to the span of h
M2 = 0; %M1-1;      % FB filter: one less than the FF filter
D = (N-1)/2;   % D is chosen large first and then decreased
K = length(x);
a_k = zeros(K,1);

% Zero padding of the i.r.
nb0 = 20;
nf0 = 20;
h = [zeros(nb0,1); h; zeros(nf0,1)];

% Get the Weiner-Hopf solution
p = zeros(M1, 1);
for i = 0:(M1 - 1)
    p(i+1) = sigma_a * conj(h(N1+nb0+1+D-i));
end

R = zeros(M1);
for row = 0:(M1-1)
    for col = 0:(M1-1)
        first_sum = (h((nb0+1):(N1+N2+nb0+1))).' * ...
            conj(h((nb0+1-(row-col)):(N1+N2+nb0+1-(row-col))));
        second_sum = (h((N1+nb0+1+1+D-col):(N1+nb0+1+M2+D-col))).' * ...
            conj((h((N1+nb0+1+1+D-row):(N1+nb0+1+M2+D-row))));
        r_w = (row == col) * sigma_w; % This is a delta only if there is no g_M.
        
        R(row+1, col+1) = sigma_a * (first_sum - second_sum) + r_w;
        
    end
end

c_opt = inv(R) * p;

b = zeros(M2,1);
for i = 1:M2
    b(i) = - (fliplr(c_opt.')*h((i+D+N1+nb0+1-M1+1):(i+D+N1+nb0+1)));
end

%% Threshold detector
x_FF = zeros(length(x),1);
for k = 0:length(x)-1
    if (k < M1 - 1)
        xconv = [flipud(x(1:k+1)); zeros(M1 - k - 1, 1)];
    else
        xconv = flipud(x(k-M1+1 + 1:k + 1));
    end
    
    x_FF(k+1) = c_opt.'*xconv;
end

detected = zeros(length(x), 1);

for i = 1:length(x_FF)
    curr = x_FF(i);
    if (real(curr) > 0)
        if (imag(curr) > 0)
            curr = 1+1i;
        else curr = 1-1i;
        end
    else
        if (imag(curr) > 0)
            curr = -1+1i;
        else curr = -1-1i;
        end
    end
    detected(i) = curr;
end

figure
for i = D + 1:length(x)
    plot(rT(i-D), 'or'), hold on, plot(x_FF(i), 'ob'), hold on,
    plot(trainingsymbols(i-D), 'ok');
    xlim([-2, 2]), ylim([-2, 2]), grid on;
    pause
    hold off
end

% TODO check that the R matrix should be Hermitian and Toeplitz for a
% LE, while it should be only Hermitian for a DFE
% Threshold to take the decision
decisions = (1+1i) * ones(K,1); % TODO set this as the actual output
%end