function  rgb_dem = demosaic_ARI2(mosaic, mask, pattern, DynamicRange)

% guided filter epsilon
eps = 1e-10;
% guidede filter window size
h = 5;
v = 5;

% green interpolation
green = green_interpolation(mosaic, mask, pattern, eps, DynamicRange);

% red and blue interpolation (first step: diagonal)
[red, blue]  =  red_blue_interpolation_first(green, mosaic, mask, eps, DynamicRange);
% red and blue interpolation (second step: horizontal/vertical)
[red, blue]  =  red_blue_interpolation_second(green, red, blue, mask, eps, DynamicRange);

% result image
rgb_dem = zeros(size(mosaic));
rgb_dem(:,:,1) = red;
rgb_dem(:,:,2) = green;
rgb_dem(:,:,3) = blue;

end


