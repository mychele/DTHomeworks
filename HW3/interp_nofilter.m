function [ out ] = interp_nofilter( x, Q )
%INTERP_NOFILTER 

out = zeros(length(x)*Q, 1);
out(1 : Q : length(out)-Q+1) = x;

end

