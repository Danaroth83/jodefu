%% Projection on the l_1,infinity norm ball of radius gamma
% 
% Description:
% This function implements the projection of x over the l_{1,inf} ball
% of radius gamma>=0 :
%   projball_l1inf(x,gamma) = argmax_{x' : l1infnorm(x')<=gamma} ||x'-x||_2
% where:
%   l1infnorm = @(x) gamma * max(sum(abs(x),1),[],2:ndims(x));
%
% Usage:
% y = projball_l1inf(x,gamma);
%
% Input:
% x: Input Matrix
% gamma: Radius of the projection ball
% 
% Output:
% y: Evaluation of the projection

function y = projball_l1inf(x, gamma)
    y=x-prox_linf1(x,gamma);
end
