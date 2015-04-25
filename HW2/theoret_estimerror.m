figure, hold on
sigma_w = 1/(4*10); %1/(T/Tc*snr);
Lvalues = [3, 7, 15, 31];
for L = Lvalues
    for N = 1:20
        n_short = mod(4-N, 4); % Num branches with a shorter filter than others
        % N_i is the number of coefficients of the filter of the i-th branch.
        N_i(1:4-n_short) = ceil(N/4);
        N_i(4-n_short + 1 : 4) = ceil(N/4) - 1;
        
        exp_deltahsqr(N) = sigma_w / (L+1) * sum(N_i .* (L+2-N_i) ./ (L+1-N_i)); %#ok<SAGROW>
    end
    plot(10*log10(exp_deltahsqr), 'DisplayName', sprintf('L=%d', L))
    legend('-DynamicLegend')
end
grid on, title('E[|\Delta h|^2]')
xlabel('N'), ylabel('E[|\Delta h|^2] [dB]')