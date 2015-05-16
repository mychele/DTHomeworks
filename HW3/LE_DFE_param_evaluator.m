%% This script performs parametric sweep on the parameters for LE and DFE, for
% different SNR

clear
close all
clc
rng default

Tc = 1;
T = 4 * Tc;
snr_vec = [6, 8, 10, 12, 14]; % dB

%% LE

printmsg_delete = ''; % Just to display progress updates

L_data = 127; % it should be ininfluent

N = 7;
M1_max = N+50;
D_max = (N-1)/2 + 10;
JminLE = zeros(length(snr_vec), M1_max, D_max);

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    
    %% Create, send and receive data with the given channel, just once for each snr
    T = 4*Tc;
    
    [packet, r_T4, ~] = txrc(L_data, snr_ch, T, Tc);
    
    % estimate the channel using the first 100 samples (4*length(ts))
    N = 7; % fixed, as specified by the teacher
    [ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);
    
    % sample to get r @ T
    init_offs = mod(m_opt, 4); % in T/4
    t0 = floor(m_opt/4); % from now consider T = 1
    T = 1;
    
    %% Detection begins
    rT = r_T4(init_offs+1:4:end); % data sampled in T
    x = rT(t0+1:t0+1+length(packet)-1)/h(N1+1).'; % data normalized by h0, starting from t0
    hi = h/h(N1+1).'; % impulse response normalized by h0
    
    N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
    for M1 = 1:M1_max;   % FF filter: equal to the span of h
        for D = 1:D_max
            % Print progress update
            printmsg = sprintf('snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            M2 = 0;      % FB filter not present in LE
            % Compute only the functional
            [JminLE(snr_i, M1, D)] = DFE_filter_jmin(packet, x, hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);
        end
    end
end

%save('jmin_LE', 'JminLE')

for i = 1:length(snr_vec)
    figure, mesh(1:D_max, 1:M1_max, reshape(abs(JminLE(i, :, :)), size(JminLE(i, :, :), 2), size(JminLE(i, :, :), 3)))
    title(strcat('Jmin for LE, snr= ', num2str(snr_vec(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin')
end

%% DFE

printmsg_delete = ''; % Just to display progress updates

N = 7;
M1_max = N+50;
D_max = M1-1;
JminDFE = zeros(length(snr_vec), M1_max, D_max);

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);
    
    %% Create, send and receive data with the given channel, just once for each snr
    T = 4*Tc;
    [packet, r_T4, ~] = txrc(L_data, snr_ch, T, Tc);
    
    % estimate the channel using the first 100 samples (4*length(ts))
    N = 7; % fixed, as specified by the teacher
    [ m_opt, h, est_sigmaw, N1, N2 ] = get_channel_info(r_T4(1:100), N, T);
    
    % sample to get r @ T
    init_offs = mod(m_opt, 4); % in T/4
    t0 = floor(m_opt/4); % from now consider T = 1
    T = 1;
    
    %% Detection begins
    rT = r_T4(init_offs+1:4:end); % data sampled in T
    x = rT(t0+1:t0+1+length(packet)-1)/h(N1+1).'; % data normalized by h0, starting from t0
    hi = h/h(N1+1).'; % impulse response normalized by h0
    
    N = N1+N2+1;    % For each symbol, we have N-1 interferers + the symbol
    
    for M1 = 1:M1_max;   % FF filter: equal to the span of h
        for D = 1:D_max
            % Print progress update
            printmsg = sprintf('for DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            
            M2 = N2 + M1 - 1 - D;  % FB filter: one less than the FF filter
            JminDFE(snr_i, M1, D) = DFE_filter_jmin(packet, x, hi, N1, N2, est_sigmaw, t0, D, M1, M2, 0);
        end
    end
end

%save('jmin_DFE', 'JminDFE')

for i = 1:length(snr_vec)
    figure, mesh(1:D_max, 1:M1_max, reshape(abs(JminDFE(i, :, :)), size(JminDFE(i, :, :), 2), size(JminDFE(i, :, :), 3)))
    title(strcat('Jmin for DFE, snr= ', num2str(snr_vec(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin')
end
