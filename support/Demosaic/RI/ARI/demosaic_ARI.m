function  rgb_dem = demosaic_ARI(mosaic, mask, pattern, DynamicRange)

% guided filter epsilon
eps = 1e-10;
% guidede filter window size
h = 5;
v = 5;

% green interpolation
green = green_interpolation(mosaic, mask, pattern, eps, DynamicRange);

% red and blue interpolation
red = red_interpolation(green, mosaic, mask, h, v, eps, DynamicRange);
blue = blue_interpolation(green, mosaic, mask, h, v, eps, DynamicRange);

% result image
rgb_dem = zeros(size(mosaic));
rgb_dem(:,:,1) = red;
rgb_dem(:,:,2) = green;
rgb_dem(:,:,3) = blue;

end


