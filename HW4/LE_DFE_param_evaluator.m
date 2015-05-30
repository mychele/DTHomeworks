%% This script performs parametric sweep on the parameters for DFE, for
% different SNR and the assigned channel

clear
close all
clc
rng default
snr_vec = [0, 5, 10, 15]; % dB

% From the assignment
t0 = 6;
N1 = 0;
N2 = 4;
N = N1 + N2 + 1;

%% DFE (parametric evaluation without channel coding)

printmsg_delete = ''; % Just to display progress updates

M1_max = N+50;
D_max = M1_max-1;
JminDFE = zeros(length(snr_vec), M1_max, D_max);
L_data = 128; % useless
bits = randi([0 1], 1, L_data);
symbols = bitmap(bits); 

for snr_i = 1:length(snr_vec)    
    snr_ch = snr_vec(snr_i);

    % Create, send and receive data with the given channel
    snrlin = 10^(snr_ch/10);
    [rcv_symb, sigma_w, h] = channel_output(symbols, snrlin);

    % Normalization!
    rcv_symb = rcv_symb(t0:end-7)/h(t0);
    hi = h(t0-N1:t0+N2)/h(t0);

    %% Receiver: compute Jmin 
    
    for M1_dfe = 1:M1_max;   % FF filter: equal to the span of h
        for D_dfe = 1:D_max
            % Print progress update
            printmsg = sprintf('for DFE, snr = %d, M1 = %d, D = %d\n', snr_ch, M1_dfe, D_dfe);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            
            M2_dfe = N2 + M1_dfe - 1 - D_dfe;  % FB filter: one less than the FF filter
            [JminDFE(snr_i, M1_dfe, D_dfe), ~] = DFE_filter(rcv_symb, hi.', N1, N2, sigma_w, D_dfe, M1_dfe, M2_dfe, false);
        end
    end
end

save('jmin_DFE', 'JminDFE')

for i = 1:length(snr_vec)
    figure, mesh(1:D_max, 1:M1_max, 10*log10(reshape(abs(JminDFE(i, :, :)), size(JminDFE(i, :, :), 2), size(JminDFE(i, :, :), 3))))
    title(strcat('Jmin for DFE, snr= ', num2str(snr_vec(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin [dB]')
end
