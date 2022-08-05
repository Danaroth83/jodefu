function mask = rectangularmask(nx,ny,nx0,ny0)
%
% Usage: out = rectangularmask(nx,ny,nx0,ny0)
%
% + nx, ny = dimensions of the output
% + nx0, ny0 = dimensions of the rectangular support
%
% Compute a rectangular binary mask in the Fourier domain.
%

mask = zeros(ny0,nx0);
omega_x = 1+mod(nx + ifftshift((0:nx0-1)-floor(nx0/2)), nx);
omega_y = 1+mod(ny + ifftshift((0:ny0-1)-floor(ny0/2)), ny);
mask(omega_y,omega_x) = 1;

end