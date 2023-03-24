%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This code is modified from the downloaded code of original guided filter paper. 
%       http://research.microsoft.com/en-us/um/people/kahe/eccv10/
%       "Guided Image Filtering", by Kaiming He, Jian Sun, and Xiaoou Tang, in TPAMI 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q] = guidedfilter_MLRI(G, R, mask, I, p, M, h, v, eps)

% image size
[hei, wid] = size(I);
% the number of the sammpled pixels in each local patch
N = boxfilter(M, h, v);
N(find(N == 0)) = 1;
% the size of each local patch; N=(2r+1)^2 except for boundary pixels.
N2 = boxfilter(ones(hei, wid), h, v);

mean_Ip = boxfilter(I.*p.*M, h, v) ./ N;
mean_II = boxfilter(I.*I.*M, h, v) ./ N;

% linear coefficients
a = mean_Ip ./ (mean_II + eps);
N3 = boxfilter(mask, h, v);
N3(find(N3 == 0)) = 1;
mean_G = boxfilter(G.*mask, h, v) ./ N3;
mean_R = boxfilter(R.*mask, h, v) ./ N3;
b = mean_R - a .* mean_G;

% weighted average
dif = boxfilter(G.*G.*mask,h,v).*a.*a + b.*b.*N3 + boxfilter(R.*R.*mask,h,v)...
    + 2*a.*b.*boxfilter(G.*mask,h,v) - 2*b.*boxfilter(R.*mask,h,v) - 2*a.*boxfilter(R.*G.*mask,h,v);
dif = dif ./ N3;
dif = dif.^0.5;
dif(find(dif<1e-3)) = 1e-3;
dif = 1./dif;
wdif = boxfilter(dif,h,v);
mean_a = boxfilter(a.*dif,h,v) ./ (wdif+1e-4);
mean_b = boxfilter(b.*dif,h,v) ./ (wdif+1e-4);

% final output
q = mean_a .* G + mean_b;

end

