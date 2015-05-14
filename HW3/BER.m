function [ pbit ] = BER( sent, detected )
% Computes the BER, it accepts symbols
%

if (length(sent) ~= length(detected))
    disp('Error in IBMAP, the sequences do not have the same length')
    return
end
sent_bit = ibmap(sent);
det_bit = ibmap(detected);

pbit = sum(abs(sent_bit - det_bit))/length(sent_bit);

end

