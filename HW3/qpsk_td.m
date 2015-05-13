function [ out ] = qpsk_td( in )
% Threshold detector for QPSK
if (real(in) > 0)
    if (imag(in) > 0)
        out = 1+1i;
    else
        out= 1-1i;
    end
else
    if (imag(in) > 0)
        out = -1+1i;
    else
        out = -1-1i;
    end
end

end

