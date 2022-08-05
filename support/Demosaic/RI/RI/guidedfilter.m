%%%%%% GuidedFilter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   - guidance image: I
%   - filtering input image: p
%   - binary mask: M
%   - local directional window radius: h, v
%   - regularization parameter: eps
%
%   This code is modified from the downloaded code of original guided filter paper. 
%       http://research.microsoft.com/en-us/um/people/kahe/eccv10/
%       "Guided Image Filtering", by Kaiming He, Jian Sun, and Xiaoou Tang, in TPAMI 2012.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q] = guidedfilter(I, p, M, h, v, eps)

% threshold parameter
th = 0.00001*255*255;

% Image size
[hei, wid] = size(I); 
% The number of the sammpled pixels in each local patch
N = boxfilter(M, h, v);
N(find(N == 0)) = 1;
% The size of each local patch; N=(2h+1)*(2v+1) except for boundary pixels.
N2 = boxfilter(ones(hei, wid), h, v);

mean_I = boxfilter(I.*M, h, v) ./ N;
mean_p = boxfilter(p, h, v) ./ N;
mean_Ip = boxfilter(I.*p, h, v) ./ N;

% The covariance of (I, p) in each local patch.
cov_Ip = mean_Ip - mean_I .* mean_p;

mean_II = boxfilter(I.*I.*M, h, v) ./ N;
var_I = mean_II - mean_I .* mean_I;
var_I(find(var_I<th)) = th;

% linear coefficients
a = cov_Ip ./ (var_I + eps);
b = mean_p - a .* mean_I;

mean_a = boxfilter(a, h, v) ./ N2;
mean_b = boxfilter(b, h, v) ./ N2;

% output
q = mean_a .* I + mean_b;

end
