function out = gaussian_filter(sigma,width,height)
% GAUSSIAN_FILTER:
%
% Usage: out = gaussian_filter(sigma,width,height)
%
% + sigma : gaussian standard deviation (frequency domain)
% + width, height : output dimensions.
%
% Compute a gaussian filter in the Fourier domain.
%
a = -floor(width/2) + (0:(width-1));
b = -floor(height/2) + (0:(height-1));
[A,B] = meshgrid(a,b);
out = ifftshift(exp(-2*(pi*sigma)^2*((A/width).^2+(B/height).^2)));
end