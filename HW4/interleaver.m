function [interleaved_bits] = interleaver(bits)
    % This function receives a sequence of bits and scrambles it
    % INPUT:
    % bits: the bits to interleave
    % OUTPUT:
    % interleaved_bits: the interleaved bits
       
    % Input should be a multiple of 14061600 = lcm(rows*columns, 64800) bits
    if (mod(14061600, length(bits)) ~= 0)
        printf('Length of the input vector should be a multiple of 14061600');
    end

    interleaved_bits = zeros(1,length(bits));
    
    rows = 31;
    columns = 35;
    
    % We work with a rowsxcolumns matrix
    for matrix = 0:(length(bits)/(rows*columns) - 1)
        curr_matrix = matrix * rows * columns;
        for col = 0:(columns-1)
            interleaved_bits(curr_matrix + col * rows + 1 : curr_matrix + col * rows + rows) = ...
                bits(curr_matrix + col + 1 : columns : curr_matrix + col + columns * rows);
        end
    end
end