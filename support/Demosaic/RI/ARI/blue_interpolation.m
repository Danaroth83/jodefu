function blue = blue_interpolation(green, mosaic, mask, h, v, eps, DynamicRange)

% sparse Laplacian filter
F = [ 0, 0,-1, 0, 0;
      0, 0, 0, 0, 0;
     -1, 0, 4, 0,-1;
      0, 0, 0, 0, 0;
      0, 0,-1, 0, 0];
% sparse Laplacian map of B and G
lap_blue = imfilter(mosaic(:,:,3), F, 'replicate');
lap_green = imfilter(green.*mask(:,:,3), F, 'replicate');

% tentative estimate generation
tentativeB = guidedfilter_MLRI(green, mosaic(:,:,3), mask(:,:,3), lap_green, lap_blue, mask(:,:,3), h, v, eps);
tentativeB = clip(tentativeB,0,DynamicRange);

% residual calculation
residualB = mask(:,:,3).*(mosaic(:,:,3)-tentativeB);

% bicubic interpolation
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
residualB = imfilter(residualB, H, 'replicate');

% add tentative estimate
blue = residualB + tentativeB;
blue = clip(blue,0,DynamicRange);

end
