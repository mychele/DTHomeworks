function plot_spectrum(signal, varargin)

    % Compute different signal analysis
    
    if (length(varargin) == 1)
        N_ar = varargin{1};
    elseif (isempty(varargin))
        N_ar = 3;
    else
        return
    end
        
    
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
    if (N_ar > 0)
        %compute variance of AR model and plot it to identify the knee
        %it's computed up to K/5 - 1
        autoc = autocorrelation(signal, N_corr);
%         useful only to know which is the knee of the sigma_w
%         upp_limit = 60;
%         sigma_w = zeros(1, upp_limit);
%         for N = 1:upp_limit
%             [~, sigma_w(N)] = arModel(N, autoc);
%         end
%         figure, plot(1:upp_limit, 10*log10(sigma_w))% Uncomment to plot sigma
        
        %compute the vector of coefficients a
        [a, sigma_w] = arModel(N_ar, autoc);
        [H, omega] = freqz(1, [1; a], K, 'whole');
    end

    clear a autoc D fir_bs_1 N N_corr S upp_limit window window_correlogram

    % Plot PSD estimate

    figure
    plot((1:K)/K, 10*log10(P_welch), 'Color', 'r', 'LineWidth', 2)
    hold on
    plot((1:K)/K, 10*log10(abs(correlogram)), 'Color', 'b', 'LineWidth', 1)
    hold on
    plot((1:K)/K, 10*log10(periodogr), 'c:')
    hold on
    if (N_ar > 0)
        plot(omega/(2*pi), 10*log10(sigma_w*abs(H).^2), 'Color', 'm', 'LineWidth', 1);
        legend('Welch', 'Correlogram', 'Periodogram', ['AR(' int2str(N_ar) ')'], 'Location', 'SouthWest')
    else
        legend('Welch', 'Correlogram', 'Periodogram', 'Location', 'SouthWest')
    end

    hold off
    axis([0, 1, -40, 10])
    title('Spectral analysis')

end