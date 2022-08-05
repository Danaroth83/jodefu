function  rgb_dem = demosaic_MLRI(mosaic, mask, pattern, sigma, DynamicRange)

% green interpolation
green = green_interpolation(mosaic, mask, pattern, sigma, DynamicRange); 
green = clip(green,0,DynamicRange);

% red and blue interpolation
red = red_interpolation(green, mosaic, mask, DynamicRange);
blue = blue_interpolation(green, mosaic, mask, DynamicRange);
red = clip(red,0,DynamicRange);
blue = clip(blue,0,DynamicRange);

% result image
rgb_dem(:,:,1) = red;
rgb_dem(:,:,2) = green;
rgb_dem(:,:,3) = blue;

end
