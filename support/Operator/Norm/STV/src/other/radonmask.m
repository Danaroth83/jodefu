function mask = radonmask(nx,ny,num_rays)
% RADONMASK compute a Radon tomography mask in the Fourier domain 
% 
% nx = mask width
% ny = mask height 
% num_rays = number of rays defining the mask support 

mask = zeros(ny,nx); 
T = linspace(0,pi,num_rays+1); 
T(end) = []; 
n = max(nx,ny); 
for theta = T
    rho = linspace(-1,1,3*n)*n; 
    x = round(rho*cos(theta)) + floor(nx/2) + 1; 
    y = round(rho*sin(theta)) + floor(ny/2) + 1; 
    id = boolean((x >= 1) .* (x <= nx) .* (y >= 1) .* (y <= ny)); 
    mask((x(id)-1)*ny+y(id)) = 1; 
end

mask = ifftshift(mask); 

end

