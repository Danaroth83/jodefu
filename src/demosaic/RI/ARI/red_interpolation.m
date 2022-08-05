function red = red_interpolation(green, mosaic, mask, h, v, eps, DynamicRange)

% sparse Laplacian filter
F = [ 0, 0,-1, 0, 0;
      0, 0, 0, 0, 0;
     -1, 0, 4, 0,-1;
      0, 0, 0, 0, 0;
      0, 0,-1, 0, 0];
% sparse Laplacian map of R and G
lap_red = imfilter(mosaic(:,:,1), F, 'replicate');
lap_green = imfilter(green.*mask(:,:,1), F, 'replicate');

% tentative estimate generation
tentativeR = guidedfilter_MLRI(green, mosaic(:,:,1), mask(:,:,1), lap_green, lap_red, mask(:,:,1), h, v, eps);
tentativeR = clip(tentativeR,0,DynamicRange);

% residual calculation
residualR = mask(:,:,1).*(mosaic(:,:,1)-tentativeR);

% bicubic interpoaltion
s_32 = s(3/2);
s_12 = s(1/2);
s_0 = s(0);
s_1 = s(1);
a = s_32^2;
b = s_32*s_12;
c = s_12^2;
d = s_0*s_12;
e = s_0*s_32;
f = s_1*s_12;
g = s_1*s_32;
H = [a, g, b, e, b, g, a;
     g, 0, f, 0, f, 0, g;
     b, f, c, d, c, f, b;
     e, 0, d, 1, d, 0, e;
     b, f, c, d, c, f, b;
     g, 0, f, 0, f, 0, g;
     a, g, b, e, b, g, a];
residualR = imfilter(residualR, H, 'replicate');

% add tentative estimate
red = residualR + tentativeR;
red = clip(red,0,DynamicRange);

end
