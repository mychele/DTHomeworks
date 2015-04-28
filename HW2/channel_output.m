function [ d, h_mean ] = channel_output( x, T, Tc, sigma_w, N_h, h_mat )
% CHANNEL_OUTPUT Generate channel output (that is the desired signal) via a
% polyphase implementation (for Nh<=4).

d = zeros(T*length(x),1);
h_used_coeff = zeros(4, length(x)); % check these dimensions
for k = 0:length(x)-1
    % Generate white noise
    w = wgn(4, 1, 10*log10(sigma_w), 'complex');
    for idx = 0:4-1 % actually idx varies between 0 and 4 since there are four branches
        if (idx < N_h)
            d(k*T + idx*Tc + 1) = h_mat(idx+1,k*T + idx*Tc+1) * x(k+1) + w(idx+1);
            % store the coefficient actually used, it will be useful later on
            h_used_coeff(idx + 1, k + 1) = h_mat(idx+1,k*T + idx*Tc+1);
        else
            d(k*T + idx*Tc + 1) = 0 + w(idx+1); % no ray, just the noise
            h_used_coeff(idx + 1, k + 1) = 0;
        end
    end
end

% No need to drop the coefficients of the transient, since the transient is zero when Nh <= 4.

% Compute the mean coefficient of impulse response of each ray over the interval of
% interest of L samples in order to compare them with the estimated impulse response.
h_mean = mean(h_used_coeff, 2);

end