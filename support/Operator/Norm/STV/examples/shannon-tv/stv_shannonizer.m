%%                    Image Shannonization
% 
% This example illustrates how the 'stvdenoise.m' function can be used to
% perform image shannonization.  
%
% Given an input image u0, we use a STVn regularization with weighted
% frequencies model to compute an image which is at the same time
% close to u0 and nicely interpolable (that is, interpolable without
% artifacts such as oscillatory patterns due to rining or aliasing). 
% 
% See details in [1], 
%
% [1] R. Abergel and L. Moisan, ``The Shannon Total Variation'', Journal of
% Mathematical Imaging and Vision, 2017.  
%
% This script makes use of data (images, convolution kernels) stored into
% the 'data' folder of this program, the following relative path must be
% adapted accordingly. 
%
relative_path_to_data = './data/'; 

%% Image Shannonization 

% load & display the input image u0 and its Fourier spectrum (see aliasing)
u0 = double(imread([relative_path_to_data,'brooklyn_bridge.tif'])); % reference image
figure('Name','input image'); imshow2(u0); 
figure('Name','spectrum of the input image (remark aliasing)'); imshow2(log(1+abs(fftshift(fft2(u0))))); 

% configure display function for the video mode
logdft = @(u)log(1+abs(fftshift(fft2(u)))); 
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; % use this function to display side-by-side the image and its spectrum in log scale (both displayed image are normlized at each iteration)

% set algorithm parameters, run algorithm, and display results
[ny,nx] = size(u0); 
n = 2; 
gam = gaussian_filter(1.4,nx,ny); 
lambda = 0.05; 
niter = 100;
verbose = true; 
video = true; 

[u,E] = shannonizer(u0,gam,n,lambda,niter,'verbose',verbose,'video',video,'displayFcn',displayFcn); 

figure('Name','Shannonized image'); imshow2(u); 
figure('Name','spectrum of the Shannonized image'); imshow2(log(1+abs(fftshift(fft2(u))))); 
figure(); plot(1:length(E),E); title('energy evolution'); xlabel('iteration'); ylabel('energy'); 

% resample the input image using spline interpolation (see spurious
% oscillations due to aliasing) 
x0 = 310; y0 = 200; width = 128; height = 128; x1 = x0+width-1; y1 = y0+height-1; z = 4;
[X,Y] = meshgrid(1+(x0:1/z:x1),1+(y0:1/z:y1));
u0_bicubic = interp2(u0,X,Y,'spline');
figure('Name','resampling of the input image (aliased)'); imshow2(u0_bicubic); 

% resample the shannonized image using spline interpolation (observe the
% absence of spurious oscillations)   
u_bicubic = interp2(u,X,Y,'spline');
figure('Name','resampling of the Shannonized image (aliased)'); imshow2(u_bicubic); 
