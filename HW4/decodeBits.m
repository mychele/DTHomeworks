function [decoded_bits] = decodeBits(bits)
    % Create the encoder
    dec = fec.ldpcdec;
    dec.DecisionType = 'Hard Decision';
    dec.OutputFormat = 'Information Part';
    dec.NumIterations = 50;
    dec.DoParityChecks = 'Yes';

    numInfoBits = dec.NumInfoBits; % Length of the info words
    
    if (mod(length(bits), numInfoBits) ~= 0)
        disp('Length of the input vector should be a multiple of 64800');
        return;
    end
    
    decoded_bits = zeros(length(bits)/2,1);
    % Iterate over the input info bits and encode them
    for idx = 0:(length(bits)/(2*numInfoBits))-1
        current_bits = bits(2*idx*numInfoBits + 1 : 2*idx*numInfoBits + 2*numInfoBits);
        decoded_bits(idx*numInfoBits+1:idx*numInfoBits + numInfoBits) = decode(dec, current_bits);
    end
end