function [ P_min, P_max, P_sorted, fft_mean] = findSine( z, window, S )

%   z1: signal for which to compute the estimate
%   window: an array containing the window to use to frame the signal
%   S: number of overlapping samples 

D = length(window);
K = length(z); % signal length
M_w = 1/D * sum(window.^2); % Power of the window
N_s = floor((K-D)/(D-S) + 1); % number of subsequences
P_per_w = zeros(K, N_s);
%figure
for s = 0:(N_s-1)
  z_s = window .* z( (s*(D-S) + 1) : (s*(D-S) + D) ); % 1.495 with index + 1
  Z_s = fft(z_s, K);
  Z_s_all(:, s + 1) = Z_s;
  P_per_w(:, s + 1) = abs(Z_s).^2/(D*M_w);
  subplot(2, 1, 1)
  plot(10*log10(abs(P_per_w(:, s+1))))
  title([int2str(s*(D-S) + 1) ' to ' int2str(s*(D-S) + D)])
  ylim([0 50])
  subplot(2, 1, 2)
  plot(real(z(s*(D-S) + 1 : s*(D-S) + D))), hold on
  plot(imag(z(s*(D-S) + 1 : s*(D-S) + D))), hold off
  pause
end
%P_welch = sum(P_per_w, 2)/N_s;

% P_sorted = for each f, take PSD(f) for each subsequence and sort them
% (independently on the subsequence: i.e. it does NOT sort spectra for
% different subsequences, indeed they're all mixed up)
P_sorted = zeros(size(P_per_w));
for i = 1:length(P_per_w)
    P_sorted(i, :) = sort(P_per_w(i, :), 2);
end

P_min = min(P_per_w, [], 2);
P_max = max(P_per_w, [], 2);

fft_mean = mean(Z_s_all, 2);

end