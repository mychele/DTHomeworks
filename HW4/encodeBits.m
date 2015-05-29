function [encoded_bits] = encodeBits(bits)
    % Create the encoder
    enc = fec.ldpcenc;
    numInfoBits = enc.NumInfoBits; % Length of the info words
    
    if (mod(length(bits), numInfoBits) ~= 0)
        disp('Length of the input vector should be a multiple of 32400');
        return;
    end
    
    encoded_bits = zeros(2*length(bits),1);
    % Iterate over the input info bits and encode them
    for idx = 0:(ceil(length(bits)/numInfoBits))-1
        current_bits = bits(idx*numInfoBits+1:idx*numInfoBits + numInfoBits);
        encoded_bits(2*idx*numInfoBits+1:2*idx*numInfoBits + 2*numInfoBits) = encode(enc, current_bits);
    end
end