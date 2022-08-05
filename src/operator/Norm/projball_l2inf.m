%% Projection on the l_2,infinity norm ball of radius gamma
% 
% Description:
% This function implements the projection of x over the l_{2,inf} ball of
% radius gamma>=0 :
%     projball_l2inf(x,gamma) = argmin_{x' : l2infnorm(x')<=gamma} ||x'-x||_2
% where
%     l2infnorm = @(x) max(sqrt(sum(x.^2,1)),[],2:ndims(x));
% Usage:
% y = projball_l2inf(x,gamma);
%
% Input:
% x: Input Matrix
% gamma: Radius of the projection ball
% 
% Output:
% y: Evaluation of the projection

function y = projball_l2inf (x, gamma)
	
    y = x ./ max(sqrt(sum(x.^2,1))/gamma,1);
    % y = x .* (gamma./max(sqrt(sum(x.^2,1)),gamma)); % only works for gamma>0
    
end