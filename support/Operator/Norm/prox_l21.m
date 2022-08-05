%% Proximal operator of the l_2,1 norm
% 
% Description:
% This function calculates the proximal operator of the l_21
% norm, defined as:
%     l2norm = @(x) gamma * sum(sqrt(sum(x.^2,1)),2:ndims(x));
% where gamma>=0
% This is equivalent to calculating y=x-projball_l2inf(x,gamma), where:
%     projball_l2inf(x,gamma) = argmin_{x' : l2infnorm(x')<=gamma} ||x_x'||_2
% is the projection of of x over the l_2,inf ball of radius gamma.
%
% Usage:
% y = prox_21(x, gamma);
%
% Input:
% x: Input Matrix, whose norm l_2 is applied on the first dimension,
%    l_1 on the the rest
% gamma: scalar multiplier for the norm
% 
% Output:
% y: Evaluation of the proximal operator

function y = prox_l21(x, gamma)
    y = x - projball_l2inf(x,gamma);
end

% %% Alternative implementation (only works for gamma>0)
% function y = prox_l21(x, gamma)
%     tmp = max(sqrt(sum(x.^2,1)),gamma);
%     y = x .* ((tmp-gamma) ./ tmp);
% end

% %% Alternative implementation (also works for gamma=0)
% function y = prox_l21(x, gamma)
% 	  tmp = sqrt(sum(x.^2,1));
% 	  tmp2 = (tmp-gamma)./tmp;
% 	  tmp2(tmp<=gamma) = 0;
% 	  y = x.* tmp2;
% end