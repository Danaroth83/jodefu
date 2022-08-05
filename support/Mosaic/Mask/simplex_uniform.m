%% SIMPLEX UNIFORM DISTRIBUTION
%
% Description:
% This function randomly generates a series of n positive whose sum is
% fixed to 1 so that it is uniformy distributed over the frontier of a n-1
% simplex
%
% Usage:
% out=simplex_uniform(n,amount);
%
% Input:
% n: Numbers of element of each set
% amount: Number of repetitions
%
% Output:
% out: A matrix (sizes: amount x n) where each row is a vector with sum
%      equal to one
%
% Reference:
% Donald B. Rubin, The Bayesian bootstrap Ann. Statist. 9, 1981, 130-134.
% http://blog.geomblog.org/2005/10/sampling-from-simplex.html

function out=simplex_uniform(n,amount)

if nargin<=1, amount=1; end

if n>0
    out=rand(amount,n-1);
    out=sort(out,2);
    out=[out,ones(amount,1)];
    out=[out(:,1),out(:,2:end)-out(:,1:end-1)];
else
    out=[];
end

end
