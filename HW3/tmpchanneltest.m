Tc = 1;
T = 4*Tc;

L = 15;
N = 10;
mlseq = MLsequence(L);

% Replace every 0 with two 0s and every 1 with two 1s to put it into the
% bitmap
mlseqdouble = zeros(1, 2*L);
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
trainingseq = [mlseqdouble, mlseqdouble(1:2*N)];
trainingsymbols = bitmap(trainingseq);

channel_output_test(trainingsymbols, T, Tc)