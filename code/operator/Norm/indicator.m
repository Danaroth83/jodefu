%%INDICATOR 
%
% Description:
% This function returns 0 if the input condition (binary) is satisfied, 
% +Infinity if not
%
% Usage:
% out=indicator(condition);

function out = indicator(condition)
    if condition==true
        out=0;
    else
        out=Inf;
    end
end

