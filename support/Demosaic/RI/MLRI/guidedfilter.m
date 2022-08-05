function [q] = guidedfilter(G, R, mask, I, p, M, h, v, eps)

% threshold parameter
th = 0.00001*255*255;

% Image size
[hei, wid] = size(I); 
% The number of the sammpled pixels in each local patch
N = boxfilter(M, h, v);
N(find(N == 0)) = 1;
% The size of each local patch; N=(2r+1)^2 except for boundary pixels.
N2 = boxfilter(ones(hei, wid), h, v);

mean_Ip = boxfilter(I.*p.*M, h, v) ./ N;
mean_II = boxfilter(I.*I.*M, h, v) ./ N;
mean_II(find(mean_II<th)) = th;

% Eqn. (5) in the paper;
a = mean_Ip ./ (mean_II + eps);

% Eqn. (6) in the paper;
N3 = boxfilter(mask, h, v);
N3(find(N3 == 0)) = 1;
mean_G = boxfilter(G.*mask, h, v) ./ N3;
mean_R = boxfilter(R.*mask, h, v) ./ N3;
b = mean_R - a .* mean_G;

mean_a = boxfilter(a, h, v) ./ N2;
mean_b = boxfilter(b, h, v) ./ N2;

% Eqn. (8) in the paper;
q = mean_a .* G + mean_b;

end
