function [ trainingsymbols ] = ts_generation( L, Nseq )
% Training sequence generation

mlseq = MLsequence(L);

% Replace every 0 with two 0s and every 1 with two 1s to put it into the
% bitmap
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

