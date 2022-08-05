function blue = blue_interpolation(green, mosaic, mask, h, v, eps)

% B interpolation
tentativeB = guidedfilter(green, mosaic(:,:,3), mask(:,:,3), h, v,eps);
residualB = ( mosaic(:,:,3) - tentativeB ) .* mask(:,:,3);
H = [1/4,1/2,1/4;1/2,1,1/2;1/4,1/2,1/4];
residualB = imfilter(residualB, H, 'replicate');
blue = residualB + tentativeB;

end




