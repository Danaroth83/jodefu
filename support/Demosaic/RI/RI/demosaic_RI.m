function  rgb_dem = demosaic_RI(mosaic, mask, pattern, sigma, DynamicRange)

% green interpolation
green = green_interpolation(mosaic, mask, pattern, sigma, DynamicRange); 

% parameters for guided upsampling
h = 5;
v = 5;
eps = 0;

% red and blue interpolation
red = red_interpolation(green, mosaic, mask, h, v, eps);
blue = blue_interpolation(green, mosaic, mask, h, v, eps);

% result image
rgb_dem(:,:,1) = red;
rgb_dem(:,:,2) = green;
rgb_dem(:,:,3) = blue;

end
