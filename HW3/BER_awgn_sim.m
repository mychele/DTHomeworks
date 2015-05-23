%% BER AWGN simulator

% Useful data
sigma_a_2 = 2;

L_data = 2^20-1;
% create a QPSK sequence
dataseq = MLsequence(L_data);
dataseq_long = [dataseq; dataseq; dataseq; dataseq; dataseq; dataseq];
datasymbols = bitmap(dataseq_long); % the ML sequence has an odd length,
% qpsk requires even symbols
%% 
snr_vec = 5:15;
numsim = 1;
pbit_AWGN_sim = zeros(length(snr_vec), numsim);
num_bit_error_AWGN_sim = zeros(length(snr_vec), numsim);
for snr_i = 1:length(snr_vec)
    snr_ch = 10^(snr_vec(snr_i)/10);
    sigma_w_2 = sigma_a_2/snr_ch;
    fprintf('snr = %d \n', 10*log10(snr_ch));
    for j = 1:numsim
        % add noise
        w = wgn(length(datasymbols), 1, 10*log10(sigma_w_2), 'complex');
        x = datasymbols + w;
        
        detected = zeros(length(x), 1);
        % perform detection
        for i = 1:length(x)
            detected(i) = qpsk_td(x(i));
        end
        
        % compute BER
        [pbit, num_bit_error] = BER(datasymbols, detected);
        pbit_AWGN_sim(snr_i, j) = pbit;
        num_bit_error_AWGN_sim(snr_i, j) = num_bit_error;
    end
end

pbit_AWGN_mean = mean(pbit_AWGN_sim, 2);

%save('BER_awgn_sim', 'pbit_AWGN_sim');

BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, pbit_AWGN_mean), hold on, semilogy(snr_vec, BER_ideal)
xlabel('snr [dB]'), ylabel('BER'), legend('AWGN, simulation', 'AWGN')
ylim([10^-8, 10^-1])