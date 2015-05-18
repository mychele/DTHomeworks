%% This script performs parametric sweep on the parameters for LE and DFE, for
% different SNR

clear
close all
clc
rng default
snr_vec = [6, 8, 10, 12, 14]; % dB

% from ex 1
assumed_m_opt = 10;
init_offs = mod(assumed_m_opt, 4);
assumed_dly = 2;
N = 5;

%% LE
printmsg_delete = ''; % Just to display progress updates
L_data = 127; % it should be ininfluent
M1_max = N+50;
D_max = (N-1)/2 + 15;
JminLE = zeros(length(snr_vec), M1_max, D_max);

for snr_i = 1:length(snr_vec)
    snr_ch = snr_vec(snr_i);

    % Create, send and receive data with the given channel
    fprintf('Generating input symbols and channel output... ')
    [packet, r, sigma_w] = txrc(L_data, snr_ch, assumed_m_opt);
    fprintf('done!\n')
    
    % Estimate the channel using the first 100 samples (4*length(ts))
    fprintf('Estimating timing phase and IR... ')
    [ h, est_sigmaw, N1, N2 ] = get_channel_info(r(init_offs+1:25+init_offs), 0, 4, assumed_m_opt);
    fprintf('done!\n')
    
    % Sample to get r @ T
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    for M1 = 1:M1_max;   % FF filter: equal to the span of h
        for D = 1:D_max
            % Print progress update
            printmsg = sprintf('snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            M2 = 0;      % FB filter not present in LE
            % Compute only the functional
            [JminLE(snr_i, M1, D)] = DFE_filter_jmin(packet, x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2, est_sigmaw, assumed_dly, D, M1, M2, 0);
        end
    end
end

save('jmin_LE', 'JminLE')

for i = 1:length(snr_vec)
    figure, mesh(1:D_max, 1:M1_max, 10*log10(reshape(abs(JminLE(i, :, :)), size(JminLE(i, :, :), 2), size(JminLE(i, :, :), 3))))
    title(strcat('Jmin for LE, snr= ', num2str(snr_vec(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin')
end

%% DFE

printmsg_delete = ''; % Just to display progress updates

M1_max = N+50;
D_max = M1-1;
JminDFE = zeros(length(snr_vec), M1_max, D_max);

for snr_i = 1:length(snr_vec)    
    snr_ch = snr_vec(snr_i);

    % Create, send and receive data with the given channel
    fprintf('Generating input symbols and channel output... ')
    [packet, r, sigma_w] = txrc(L_data, snr_ch, assumed_m_opt);
    fprintf('done!\n')
    
    % Estimate the channel using the first 100 samples (4*length(ts))
    fprintf('Estimating timing phase and IR... ')
    [ h, est_sigmaw, N1, N2 ] = get_channel_info(r(init_offs+1:25+init_offs), 0, 4, assumed_m_opt);
    fprintf('done!\n')
    
    % Sample to get r @ T
    x = r / h(N1+1).';         % data normalized by h0
    hi = h / h(N1+1).';         % impulse response normalized by h0
    
    for M1 = 1:M1_max;   % FF filter: equal to the span of h
        for D = 1:D_max
            % Print progress update
            printmsg = sprintf('for DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1, D);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            
            M2 = N2 + M1 - 1 - D;  % FB filter: one less than the FF filter
            JminDFE(snr_i, M1, D) = DFE_filter_jmin(packet, x(1+assumed_dly : assumed_dly+length(packet)), hi, N1, N2, est_sigmaw, assumed_dly, D, M1, M2, 0);
        end
    end
end

save('jmin_DFE', 'JminDFE')

for i = 1:length(snr_vec)
    figure, mesh(1:D_max, 1:M1_max, 10*log10(reshape(abs(JminDFE(i, :, :)), size(JminDFE(i, :, :), 2), size(JminDFE(i, :, :), 3))))
    title(strcat('Jmin for DFE, snr= ', num2str(snr_vec(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin')
end
