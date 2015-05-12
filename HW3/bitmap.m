function [output] = bitmap(input)
    % Check if the input array has even length
    L = length(input);
    if (mod(L, 2) ~= 0)
        disp('Must input an even length array');
        return;
    end
    
    output = zeros(L,1);
    
    % Map each couple of values to the corresponding symbol
    for idx = 1:2:L-1
        if (isequal(input(idx:idx+1), [0; 0] ))
            output(idx) = -1-1i;
        elseif (isequal(input(idx:idx+1), [0; 1] ))
            output(idx) = 1-1i;
        elseif (isequal(input(idx:idx+1), [1; 0] ))
            output(idx) = -1+1i;
        elseif (isequal(input(idx:idx+1), [1; 1] ))
            output(idx) = +1+1i;
        end
    end
    
    % Finally, only keep the useful values of the output
    output = output(1:2:end);
    
end