%% Proximity operator of the l_1 norm
% 
% Description:
% This function calculates the proximal operator of the l_1
% norm, defined as:
%   l1norm = @(x) gamma * sum(abs(x(:))) where gamma>=0
% This is equivalent to calculating y=x-projball_linf(x,gamma), where:
%   projball_linf(x,gamma) = argmin_{x' : linf_norm(x')<=gamma} ||x_x'||_2
% is the projection of of x over an l_inf ball of radius gamma
%
% Usage:
% y = prox_l1(x, gamma);
%
% Input:
% x: Input Matrix, whose norm l_1 is applied on the first dimension,
%    l_infinity on the the second
% gamma: scalar multiplier for the norm
% 
% Output:
% y: Evaluation of the proximal operator


function y = prox_l1(x, gamma)

    y = sign(x).*max(abs(x)-gamma,0);
    % y = x - projball_linf(x,gamma); % slower
    
end