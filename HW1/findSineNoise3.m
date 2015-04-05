function [ freq ] = findSineNoise3( x, resolutionMultiplier )


% Don't use the moving average
useMA = false;
MOV_AVG_LEN = 6;

% Window
%x = x .* triang(length(x));
x = x .* kaiser(length(x), 5.65);

% Base resolution is such that a main lobe is defined by 5 frequency
% samples (4 intervals)
mainLobeWidthSamples = 4;
baseResolution = 0.5 * length(x) * mainLobeWidthSamples;
baseResolution = round(baseResolution / 2); % half the resolution, since with triang windows the lobe's width is doubled
resolution = round(baseResolution * resolutionMultiplier);
mainLobeWidthSamples = mainLobeWidthSamples * resolutionMultiplier;

% Compute FFT with a certain resolution (using zero-padding)
paddingLength = resolution - length(x);
X = fft([x; zeros(paddingLength, 1)]);
f = (0 : resolution-1) / resolution;

if (useMA)
    % Define a moving average to apply to the DFT
    movingAvgSize = MOV_AVG_LEN * mainLobeWidthSamples;
    delay = round(movingAvgSize / 2 + 1);
    mavgOfFreqSpectrum = filter(ones(1, movingAvgSize)/movingAvgSize, 1, repmat(abs(X), 3, 1));
    mavgOFSlen = length(mavgOfFreqSpectrum) / 3;
    mavgOfFreqSpectrum = mavgOfFreqSpectrum(mavgOFSlen+delay : 2*mavgOFSlen+delay-1)';
    
    % Compute the distance between the DFT and its moving average
    Y = abs(abs(X)' - mavgOfFreqSpectrum);
else
    Y = abs(X)';
end

% Plots
%figure, plot(f, abs(X))
%hold on, plot(f, mavgOfFreqSpectrum, 'r')
%plot(f, Y, 'g')


% Now find the peaks in the remaining part.
[ppks, plocs] = findpeaks(Y);
[npks, nlocs] = findpeaks(-Y);
npks = -npks;



locs = [plocs', 10*log10(ppks'),  ones(size(plocs))'; ...
        nlocs', 10*log10(npks'), -ones(size(nlocs))'];

THRESH = 1;
locs = sortrows(locs);
for i = 2:size(locs, 1)
    if (locs(i, 3) * locs(i-1, 3) == 0)
        continue;
    end
    if (abs(locs(i, 2) - locs(i-1, 2)) < THRESH)
        if (locs(i, 3) + locs(i-1, 3) == 0)
            locs(i, 3)   = 0;
            locs(i-1, 3) = 0;
        else
            locs(i, 3) = 0;
            locs(i-1, 1) = round(mean(locs((i-1):1, 1)));
        end
    end
    
end
locs = locs(find(locs(:, 3)), :);
    

plocs = locs(find(locs(:, 3) > 0), 1);
nlocs = locs(find(locs(:, 3) < 0), 1);
ppks = locs(find(locs(:, 3) > 0), 2);
npks = locs(find(locs(:, 3) < 0), 2);

% UNCOMMENT THIS  TO SEE WHAT HAPPENS
% hold off, plot(f, 10*log10(Y), 'r')
% hold on
% stem(plocs/resolution, ppks)
% stem(nlocs/resolution, npks, 'x')

freq = plocs/resolution;

end