function [p,s] = perdecomp(u)

%% Periodic plus Smooth Image Decomposition
%
% author: Lionel Moisan (Scilab language)
%       
% This function computes the periodic (p) and smooth (s) components
% of an image (2D array) u
% usage:    p = perdecomp(u)    or    [p,s] = perdecomp(u)
%
% note: this function also works for 1D signals (line or column vectors)
%
% v1.0 (07/2012): initial Scilab version
% v1.1 (10/2012): simplified coordinates, use of meshgrid
% v1.2 (10/2012): removed slow meshgrid
% v2.0 (12/2016): adaptation to Matlab language (Rémy Abergel)
  
[ny,nx,nz] = size(u);
X = 1:nx; Y = 1:ny;
v = zeros(size(u));
v(1,X,:)  = u(1,X,:)-u(ny,X,:);
v(ny,X,:) = -v(1,X,:);
v(Y,1,:) = v(Y,1,:)+u(Y,1,:)-u(Y,nx,:);
v(Y,nx,:) = v(Y,nx,:)-u(Y,1,:)+u(Y,nx,:);
fx = repmat(ones(ny,1)*cos(2.*pi*(X-1)/nx),[1,1,nz]);
fy = repmat(cos(2.*pi*(Y'-1)/ny)*ones(1,nx),[1,1,nz]);
fx(1,1,:)=0.; % avoid division by 0 in the line below
s = real(ifft2(fft2(v)*0.5./(2.-fx-fy)));
p = u-s;
end
