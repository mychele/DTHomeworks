function plot_spectrum(signal)

    % Compute different signall analysis
    
    K = length(signal);

    %PERIODOGRAM pg 84
    Z = fft(signal);
    periodogr = abs(Z).^2/K;

    %compute WELCH estimator pg 85 
    D = 200; % window size
    window = kaiser(D, 5.65);
    S = D/2; %common samples
    P_welch = welchPsd(signal, window, S);

    %CORRELOGRAM
    N_corr = ceil(K/5); % N_corr is the order of the autocorrelation estimate
    window_correlogram = kaiser(2*N_corr + 1, 5.65); % window centered around N_corr
    correlogram = correlogramPsd(signal, window_correlogram, N_corr);

    %AR model: order estimation
    %compute variance of AR model and plot it to identify the knee
    %it's computed up to K/5 - 1, don't know if it makes sense
    autoc = autocorrelation(signal, N_corr);
    upp_limit = 60;
    sigma_w = zeros(1, upp_limit);
    for N = 1:upp_limit
        [~, sigma_w(N)] = arModel(N, autoc);
    end
    %figure, plot(1:upp_limit, 10*log10(sigma_w))% Uncomment to plot sigma

    %the knee is apparently at N = 3
    %compute the vector of coefficients a
    N = 3;
    [a, sigma_w] = arModel(N, autoc);
    [H, omega] = freqz(1, [1; a], K, 'whole');

    clear a autoc D fir_bs_1 N N_corr S upp_limit window window_correlogram

    % Plot PSD estimate

    figure
    plot((1:K)/K, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
    hold on
    plot((1:K)/K, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
    hold on
    plot((1:K)/K, 10*log10(periodogr), 'c:')
    hold on
    plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
    hold off
    axis([0, 1, -40, 10])
    legend('Welch', 'Correlogram', 'Periodogram', 'AR(3)', 'Location', 'SouthEast')
    title('Spectral analysis')

end