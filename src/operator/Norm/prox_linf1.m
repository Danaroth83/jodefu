%% Proximity operator of the l_infinity,1 norm
% 
% Description:
% This function calculates the proximal operator of the l_{infinity,1}
% norm, defined as:
%   linf1norm = @(x) gamma * sum(max(abs(x),[],1),2:ndims(x)), where gamma>=0
% This is equivalent to calculating y=x-projball_l1inf(x,gamma), where:
%   projball_l1inf(x,gamma) = argmin_{x' : l1inf_norm(x')<=gamma} ||x_x'||_2
% is the projection of of x over the l_{1,inf} ball of radius gamma
%
% Usage:
% y = prox_linf1(x, gamma);
%
% Input:
% x: Input Matrix, whose norm l_1 is applied on the first dimension,
%    l_infinity on the the second, l_1 on the rest
% gamma: scalar multiplier for the norm
% 
% Output:
% y: Evaluation of the proximal operator

function y = prox_linf1(x, gamma)
    lambda = max(max((cumsum(sort(abs(x),1,'descend'),1)-gamma) ./ (1:size(x,1))',[],1),0);
    y = max(min(x, lambda), -lambda);
end

%% Alternative Implementation (slower)

% function y = prox_linf1(x, gamma)
% 
%   % Preparing the input (Allows to deal with l_inf,1 norms)
%     
%   sizex=size(x);
%   x=reshape(x,size(x,1),[]);
%     
%   % Computation
%     
% 	y = zeros(size(x));
% 	mask = sum(abs(x),1)>gamma;
% 	x = x(:,mask);
% 	lambda = max((cumsum(sort(abs(x),1,'descend'),1)-gamma) ./ (1:size(x,1))',[],1);
% 	y(:,mask) = max (min(x, lambda),-lambda);
%     
%   % Reshaping the matrix
%   y=reshape(y,sizex);
%     
% end

