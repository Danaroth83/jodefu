%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This code is modified from the downloaded code of original guided filter paper. 
%       http://research.microsoft.com/en-us/um/people/kahe/eccv10/
%       "Guided Image Filtering", by Kaiming He, Jian Sun, and Xiaoou Tang, in TPAMI 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q] = guidedfilter(I, p, M, h, v, eps)

% image size
[hei, wid] = size(I); 
% the number of the sammpled pixels in each local patch
N = boxfilter(M, h, v);
N(find(N == 0)) = 1;
% the size of each local patch; N=(2r+1)^2 except for boundary pixels.
N2 = boxfilter(ones(hei, wid), h, v);

mean_I = boxfilter(I.*M, h, v) ./ N;
mean_p = boxfilter(p.*M, h, v) ./ N;
mean_Ip = boxfilter(I.*p.*M, h, v) ./ N;

cov_Ip = mean_Ip - mean_I .* mean_p;
mean_II = boxfilter(I.*I.*M, h, v) ./ N;
var_I = mean_II - mean_I .* mean_I;

% linear coefficients
a = cov_Ip ./ (var_I + eps);
b = mean_p - a .* mean_I;

% weighted average 
dif = boxfilter(I.*I.*M,h,v).*a.*a + b.*b.*N + boxfilter(p.*p.*M,h,v)...
    + 2*a.*b.*boxfilter(I.*M,h,v) - 2*b.*boxfilter(p.*M,h,v) - 2*a.*boxfilter(p.*I.*M,h,v);
dif = dif ./ N;
dif = dif.^0.5;
dif(find(dif<1e-3)) = 1e-3;
dif = 1./dif;
wdif = boxfilter(dif,h,v);
mean_a = boxfilter(a.*dif,h,v) ./ (wdif+1e-4);
mean_b = boxfilter(b.*dif,h,v) ./ (wdif+1e-4);

% final output
q = mean_a .* I + mean_b;

end
