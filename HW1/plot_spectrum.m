function plot_spectrum(signal, N_ar)
% This function computes different estimators on a given signal: Periodogram, 
% Welch periodogram, Correlogram and if N_ar > 0 also the AR model of order
% N_ar. The parameters and the windows of Welch and Correlogram are hard
% coded in the function and not passed as arguments. It also plots the
% various estimate on the same plot, with frequency normalized in [0,1].
% Use ylim in the main script to set the desired dinamic in the Y axis.

    K = length(signal);

    % PERIODOGRAM pg 84
    Z = fft(signal);
    periodogr = abs(Z).^2/K;

    % compute WELCH estimator pg 85 
    D = 200; % window size
    window = kaiser(D, 5.65);
    S = D/2; %common samples
    P_welch = welchPsd(signal, window, S);

    % CORRELOGRAM
    N_corr = ceil(K/5); % N_corr is the order of the autocorrelation estimate
    window_correlogram = kaiser(2*N_corr + 1, 5.65); % window centered around N_corr
    correlogram = correlogramPsd(signal, window_correlogram, N_corr);

    % AR model
    if (N_ar > 0)
        %compute variance of AR model and plot it to identify the knee
        %it's computed up to K/5 - 1
        autoc = autocorrelation_biased(signal, N_corr);
        
        %compute the vector of coefficients a
        [a, sigma_w] = arModel(N_ar, autoc);
        [H, omega] = freqz(1, [1; a], K, 'whole');
    end

    
    % Plot PSD estimate
    figure, hold on
    plot((1:K)/K, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
    plot((1:K)/K, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
    plot((1:K)/K, 10*log10(periodogr), 'c:')
    if (N_ar > 0)
        plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
        legend('Welch', 'Correlogram', 'Periodogram', ['AR(' int2str(N_ar) ')'], 'Location', 'SouthWest')
    else
        legend('Welch', 'Correlogram', 'Periodogram', 'Location', 'SouthWest')
    end
    hold off
    title('Spectral analysis')
    xlabel('Normalized frequency')
    ylabel('Magnitude (dB)')

end