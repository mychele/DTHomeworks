function [output] = ibmap(input)
    % Check if the input array has even length
    L = length(input);
    
    output = zeros(2*L,1);
    
    % Map each couple of values to the corresponding symbol
    % The real part gives the bit 
    for k = 1:2:length(output)-1
        symbol = input((k+1)/2);
        if (real(symbol) == 1)
           b2k = 1;
        else
           b2k = 0;
        end
        
        if (imag(symbol) == 1)
            b2k1 = 1;
        else
            b2k1 = 0;
        end
        
        output(k)= b2k;
        output(k+1)= b2k1;
    end    
end