function  rgb_dem = demosaic_MLRI2(mosaic, mask, pattern, sigma, eps, DynamicRange)

% green interpolation
green = green_interpolation(mosaic, mask, pattern, sigma, eps, DynamicRange); 
green = clip(green,0,DynamicRange);

% red and blue demosaicking
red = red_interpolation(green, mosaic, mask, eps, DynamicRange);
blue = blue_interpolation(green, mosaic, mask, eps, DynamicRange);
red = clip(red,0,DynamicRange);
blue = clip(blue,0,DynamicRange);

% result image
rgb_dem(:,:,1) = red;
rgb_dem(:,:,2) = green;
rgb_dem(:,:,3) = blue;

end
