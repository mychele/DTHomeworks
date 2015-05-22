%% BER AWGN simulator

% Useful data
q = [0,0,0,0,0,0,0,0,0.19*exp(-1i*2.21), 0.09*exp(1i*1.64), 0.7*exp(-1i*2.57), ...
    0.45, 0.6*exp(-1i*2.26), 0.35*exp(1i*3.15), 0.24*exp(1i*1.34), 0.37*exp(1i*2.6), ...
    0.34*exp(-1i*1.17), 0, 0.15*exp(-1i*2.66), 0.15*exp(1i*3.27), 0.17*exp(1i*2.13), ...
    0.4*exp(1i*2.06), 0.58*exp(-1i*1.51), 0.03*exp(1i*2.15), 0.18*exp(1i*3.6), ...
    0.29*exp(1i*3.17), 0.4*exp(-1i*1.63), 0.07*exp(-1i*3.16)];
h = q(3:4:end);
E_h = sum(abs(h).^2);
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
    sigma_w_2 = sigma_a_2*E_h/snr_ch;
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


BER_ideal = BER_awgn(snr_vec);

figure, semilogy(snr_vec, pbit_AWGN_mean), hold on, semilogy(snr_vec, BER_ideal)
xlabel('snr [dB]'), ylabel('BER'), legend('AWGN, simulation', 'AWGN')
ylim([10^-8, 10^-1])