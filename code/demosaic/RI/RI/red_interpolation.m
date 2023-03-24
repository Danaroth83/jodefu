function red = red_interpolation(green, mosaic, mask, h, v, eps)

% R interpolation
tentativeR = guidedfilter(green, mosaic(:,:,1), mask(:,:,1), h, v, eps);
residualR = ( mosaic(:,:,1) - tentativeR ) .* mask(:,:,1);
H = [1/4,1/2,1/4;1/2,1,1/2;1/4,1/2,1/4];
residualR = imfilter(residualR, H, 'replicate');
red = residualR + tentativeR;

end




