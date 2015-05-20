function [ p ] = MLsequence( L )
% Generate a Maximum Length Pseudo Noise sequence, using shift and xor
% operators. L is the desired length of the resulting PN sequence.

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
        case 2^10 -1
            p(l) = xor(p(l-7), p(l-10));
        case 2^11 -1
            p(l) = xor(p(l-9), p(l-11));
        case 2^15 -1
            p(l) = xor(p(l-14), p(l-15));
        case 2^18 -1
            p(l) = xor(p(l-11), p(l-18));
        case 2^19 -1
            p(l) = xor(xor(xor(p(l-14), p(l-17)), p(l-18)), p(l-19));
        case 2^20 -1
            p(l) = xor(p(l-17), p(l-20));
        otherwise
            p = [];
            disp('MLsequence: length not supported');
    end
end
end