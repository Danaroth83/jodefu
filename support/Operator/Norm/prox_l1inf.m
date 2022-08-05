%% Proximity operator of the l_1,infinity norm
% 
% Description:
% This function implements the algorithm described in [1] for calculating
% the proximal operator of the l_1,infinity norm, defined as:
%   l1inf1norm = @(x) gamma * max(sum(abs(x),1),[],2:ndims(x)); 
% with gamma>=0
% This is equivalent to calculating y=x-projball_linf1(x,gamma), where:
%   projball_linf1inf(x,gamma) = argmin_{x' : linf1infnorm(x')<=gamma} ||x_x'||_2
% is the projection of of x over the l_inf,1,inf ball of radius gamma).
%
% Usage:
% y = prox_l1inf(x, gamma);
%
% Input:
% x: Input Matrix, whose norm l_1 is applied on the first dimension,
%    l_inf on the second, l_1 on the rest
% gamma: scalar multiplier for the norm
% 
% Output:
% y: Evaluation of the proximal operator
%
% References:
% [1] L. Condat, "Fast projection onto the simplex and the l1 ball",
% Mathematical Programming, vol. 158, no. 1, pp. 575-585, 2016. 

function y = prox_l1inf(x, gamma)
    y=x-projball_linf1(x,gamma);
end
