%% Projection on the l_2,1,infinity norm ball of radius gamma
% 
% Description:
% This function implements the projection of x over the l_{2,1,inf} ball of
% radius gamma>=0 :
%     projball_l21inf(x,gamma) = argmin_{x' : l2infnorm(x')<=gamma} ||x'-x||_2
% where
%     l21infnorm = @(x) max(sum(sqrt(sum(x.^2,1)),2),[],3:ndims(x));
% Usage:
% y = projball_l2inf(x,gamma);
%
% Input:
% x: Input Matrix
% gamma: Radius of the projection ball
% 
% Output:
% y: Evaluation of the projection

function y = projball_l21inf(x, gamma)
    y=x-prox_l2inf1(x,gamma);
end