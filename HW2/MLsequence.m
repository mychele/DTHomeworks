function [ p ] = MLsequence( L )
%MLSEQUENCE

r = log2(L+1);
p = zeros(L,1);
p(1:r) = ones(1,r).'; % Set arbitrary initial condition
for l = r+1:(L)
    switch L
        case 3
            p(l) = xor(p(l-1), p(l-2));
        case 7
            p(l) = xor(p(l-2), p(l-3));
        case 15
            p(l) = xor(p(l-3), p(l-4));
        case 31
            p(l) = xor(p(l-3), p(l-5));
        case 63
            p(l) = xor(p(l-5), p(l-6));
        case 127
            p(l) = xor(p(l-6), p(l-7));
    end
end
end