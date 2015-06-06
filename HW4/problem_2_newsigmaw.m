% This script solves the second problem.
clear, clc, close all
rng default

% Initialize parameters based on the assigned channel
t0 = 5;
N1 = 0;
N2 = 4;
M1_dfe = 15;
D_dfe = M1_dfe - 1;
M2_dfe = N2 + M1_dfe - 1 - D_dfe;
sigma_a = 2;

parpool(15);

%% Estimated channel, DFE, coded data
% In this section we don't use the txrc function for now. This is because
% txrc only accepts L_data such that an ML sequence can be directly
% created.

L = 31;
Nseq = 7;

% Get optimal number of bits
desired_bits = 2^24;
% Compute the closest number of bits that both interleaver and encoder will like
search_step = 32400;
bit_number = ceil(desired_bits / search_step) * search_step;

numsim = 10;

snr_vec_estch_coded_new = [1, 2, 3:0.02:3.6];  % Pbit falls at 3.5 dB
seq_lengths_estch_coded_new = bit_number*ones(1, length(snr_vec_estch_coded_new));
Pbit_estch_coded_new = zeros(length(snr_vec_estch_coded_new),numsim);

for sim = 1:numsim
    parfor snr_idx = 1:length(snr_vec_estch_coded_new)
        curr_snr = snr_vec_estch_coded_new(snr_idx);
        fprintf('Estimated channel, coded, snr = %.2f\n', curr_snr);
        
        % Generate the current needed sequence
        packet = randi([0 1], 1, seq_lengths_estch_coded_new(snr_idx));
        
        enc_packet = encodeBits(packet);
        int_enc_packet = interleaver(enc_packet);  % Interleave the encoded bits
        
        symbols = [ts_generation(L, Nseq); bitmap(int_enc_packet.')];
        
        % Send through the channel
        [rcv_symb, sigma_w, ~] = channel_output(symbols, 10^(curr_snr/10), false);
        
        % Perform estimation
        [h, est_sigma_w] = get_channel_info(rcv_symb(t0+1:t0+L+Nseq), N1, N2);
        
        rcv_symb = rcv_symb(t0+1:end-7)/h(N1+1);
        hi = h / h(N1+1);
        
        % Receiver: filter with DFE
        [Jmin, psi, rcv_symb] = DFE_filter(rcv_symb, hi, N1, N2, est_sigma_w, D_dfe, M1_dfe, M2_dfe, true, false);
        
        noise_var = (Jmin-sigma_a*abs(1-psi(D_dfe+1))^2)/abs(psi(D_dfe+1))^2; % This includes white noise and isi
        
        % Compute Log Likelihood Ratio
        llr = zeros(2*length(packet),1);
        llr(1:2:end) = -2*real(rcv_symb(L+Nseq+1:end))/(noise_var/2);
        llr(2:2:end) = -2*imag(rcv_symb(L+Nseq+1:end))/(noise_var/2);
        
        llr = deinterleaver(llr); % Deinterleave the loglikelihood ratio first
        
        dec_packet = decodeBits(llr).';
        
        % Compute the Pbit and store it
        Pbit_estch_coded_new(snr_idx, sim) = sum(xor(dec_packet, packet))/length(packet);
        
    end
end
% Save current results
save('Problem2_estch_coded_new', 'snr_vec_estch_coded_new', ...
    'seq_lengths_estch_coded_new', 'Pbit_estch_coded_new');

%% Plot BER graphs

% load('Problem2_estch_uncoded.mat');
% load('Problem2_estch_coded.mat');
% load('Problem2_knownch_uncoded.mat');
% load('Problem2_knownch_coded.mat');
% 
% figure, semilogy(snr_vec_knownch_uncoded, Pbit_knownch_uncoded), hold on
% semilogy(snr_vec_knownch_coded, Pbit_knownch_coded), hold on
% semilogy(snr_vec_estch_uncoded, Pbit_estch_uncoded, '--')
% semilogy(snr_vec_estch_coded, Pbit_estch_coded, '--')
% xlabel('\Gamma [dB]'), ylabel('Pbit'), grid on
% ylim([10^-5, 10^-1])
% xlim([0, 14])
% legend('Known channel, uncoded', 'Known channel, coded', ...
%     'Estimated channel, uncoded', 'Estimated channel, coded');
% title('Bit Error Rate for a DFE receiver')

%% Clean parpool
delete(gcp)