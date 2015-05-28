function [ trainingsymbols ] = ts_generation( L, Nseq )
% Training sequence generation
% This function outputs a partially repeated ML sequence in which the 
% symbols are 1+j and -1-j (two orthogonal symbols from the QPSK alphabet).

mlseq = MLsequence(L); % Get the 0-1 ML sequence

% Replace every 0 with two 0s and every 1 with two 1s 
mlseqdouble = zeros(2*L,1);
for i = 1:L
    switch mlseq(i)
        case 0
            mlseqdouble(2*i-1) = 0;
            mlseqdouble(2*i) = 0;
        case 1
            mlseqdouble(2*i-1) = 1;
            mlseqdouble(2*i) = 1;
    end
end

% Repeat the sequence and bitmap it to get the symbols
trainingseq = [mlseqdouble; mlseqdouble(1:2*Nseq)];
trainingsymbols = bitmap(trainingseq);

end