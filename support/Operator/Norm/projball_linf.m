%% Projection on the l_infinity norm ball of radius gamma
% 
% Description:
% This function implements the projection of x over the l_inf ball of
% radius gamma>=0 :
%   projball_linf(x,gamma) = argmin_{x' : linfnorm(x')<=gamma} ||x'-x||_2
% where:
%   linfnorm= @(x) max(abs(x(:)));
%
% Usage:
% y = projball_linf(x,gamma);
%
% Input:
% x: Input Matrix
% gamma: Radius of the projection ball
% 
% Output:
% y: Evaluation of the projection

function y=projball_linf(x,gamma)

    y=max(min(x, gamma),-gamma);
    % y=sign(x).*min(abs(x),gamma);      % slower
    % y = x./max(abs(x)/gamma,1);        % slower (but still works with gamma=0)
    % y = x.*(gamma./max(abs(x),gamma)); % only works if gamma>0

end