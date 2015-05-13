function [ p ] = MLsequence( L )
% Generate a Maximum Length Pseudo Noise sequence, using shift and xor
% operators. L is the desired length of the resulting PN sequence. The
% script can handle L = 3, 7, 15, 31, 63 and 127.

% Extract r from the given L
r = log2(L+1);
p = zeros(L,1);
p(1:r) = ones(r,1); % Set arbitrary initial condition
for l = r+1:(L) % Skip the initial condition for the cycle
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
        case 2^20 -1
            p(l) = xor(p(l-17), p(l-20));
    end
end
end