%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This code is modified from the downloaded code of original guided filter paper. 
%       http://research.microsoft.com/en-us/um/people/kahe/eccv10/
%       "Guided Image Filtering", by Kaiming He, Jian Sun, and Xiaoou Tang, in TPAMI 2013.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q] = guidedfilter_diagonal(I, p, M, h, v, eps)

% diagonal window setting
r = h+v;
F = ones(2*r+1,2*r+1);
w = 2*r+1;
for i=1:v
	for t=1:2*i-1
		F(t,2*i-t) = 0;
		F(w+1-t,w+1-2*i+t) = 0;
	end
end
for i=1:h
	for t=1:2*i-1
		F(t,w+1-2*i+t) = 0;
		F(w+1-t,2*i-t) = 0;
	end
end
F2 = zeros(2*r+1,2*r+1);
F2(1:2:end,1:2:end) = 1;
F2(2:2:end,2:2:end) = 1;
F = F.*F2;

% image size
[hei, wid] = size(I); 
% the number of the sammpled pixels in each local patch
N = imfilter(M, F, 'replicate');
N(find(N == 0)) = 1;
% the size of each local patch; N=(2r+1)^2 except for boundary pixels.
N2 = imfilter(ones(hei, wid), F, 'replicate');

mean_I = imfilter(I.*M, F, 'replicate') ./ N;
mean_p = imfilter(p.*M, F, 'replicate') ./ N;
mean_Ip = imfilter(I.*p.*M, F, 'replicate') ./ N;

cov_Ip = mean_Ip - mean_I .* mean_p;
mean_II = imfilter(I.*I.*M, F, 'replicate') ./ N;
var_I = mean_II - mean_I .* mean_I;

% linear coefficients
a = cov_Ip ./ (var_I + eps);
b = mean_p - a .* mean_I;

% weighted average 
dif = imfilter(I.*I.*M,F,'replicate').*a.*a + b.*b.*N + imfilter(p.*p.*M,F,'replicate')...
    + 2*a.*b.*imfilter(I.*M,F,'replicate') - 2*b.*imfilter(p.*M,F,'replicate')...
    - 2*a.*imfilter(p.*I.*M,F,'replicate');
dif = dif ./ N;
dif = dif.^0.5;
dif(find(dif<1e-3)) = 1e-3;
dif = 1./dif;
wdif = imfilter(dif,F,'replicate');
mean_a = imfilter(a.*dif,F,'replicate') ./ (wdif+1e-4);
mean_b = imfilter(b.*dif,F,'replicate') ./ (wdif+1e-4);

% final output
q = mean_a .* I + mean_b;

end

