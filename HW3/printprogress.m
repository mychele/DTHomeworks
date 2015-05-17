function [ msg_delete ] = printprogress( msg, msg_delete )
%PRINTPROGRESS Summary of this function goes here
%   Detailed explanation goes here

% Print progress update
fprintf([msg_delete, msg]);
msg_delete = repmat(sprintf('\b'), 1, length(msg));

end

