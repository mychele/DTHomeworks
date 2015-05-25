function [ pbit, num_bit_error ] = BER( sent, detected )
% Computes the BER, it accepts symbols

if (length(sent) ~= length(detected))
    disp('Error in IBMAP, the sequences do not have the same length')
    return
end
sent_bit = ibmap(sent);
det_bit = ibmap(detected);

num_bit_error = sum(abs(sent_bit - det_bit));

pbit = num_bit_error/length(sent_bit);

end

