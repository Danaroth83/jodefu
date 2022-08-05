%%          Discrete Total Variation based image deconvolution
% 
% This example illustrates how the generic 'solver_tv_relax.m' function
% can be used to perform image deconvolution. 
%
% Given an input image u0, compute a minimizer of the energy
%
%           E(u) := ||A(u)-u0||^2 +lambda*TV(u)
%         
% where A(u) denotes the convolution between u and a given kernel (see [1]
% for the exact definition)
%
% [1] R. Abergel and L. Moisan, ``The Shannon Total Variation'', Journal of
% Mathematical Imaging and Vision, 2017. 
%
% This script makes use of data (images, convolution kernels) stored into
% the 'data' folder of this program, the following relative path must be
% adapted accordingly. 
%
relative_path_to_data = './data/'; 

%% Image deconvolution

% define operator A, its adjoint adj_A, and an upper bound L_A for |||A|||
ker = imread([relative_path_to_data,'motion_kernel_1.tif']); % motion blur kernel 
A = @(u)conv2(u,ker,'valid'); % convolution operator
adj_A = @(v)conv2(v,ker(end:-1:1,end:-1:1),'full'); % adjoint of A (warning: this is only valid when ker has odd x odd dimensions)
L_A = sum(abs(ker(:))); % upper bound for the l2 induced norm of A

% load & display the reference image and its blurry and noisy version
ref = double(imread([relative_path_to_data,'clock2.tif'])); ref = ref(51:360,1:310); % reference image 
u0 = A(ref)+2*randn(size(A(ref))); % blurry and noisy image
figure('Name','reference image'); imshow2(ref); 
figure('Name','blurry and noisy'); imshow2(u0); 

% set parameters, in particular compute init = u0 extended to the full
% domain (that of ref) with pixel recopy boundary condition (this
% dramatically improves the convergence speed, compared to the default
% setting of 'init')
[ny0,nx0] = size(u0);
[ky,kx] = size(ker); 
dx = floor(kx/2); 
dy = floor(ky/2); 
init = u0([ones(1,dy),1:ny0,ny0*ones(1,dy)],[ones(1,dx),1:nx0,nx0*ones(1,dx)]);
lambda = 0.4; 
niter = 200;
verbose = true; 
video = true; 
crop_roi = @(u)u((1+dy):(end-dy),(1+dx):(end-dx)); % display function for the video mode: crop the image to the domain of u0 

% Run algorithm and display result. Notice that since the restored image
% has a bigger domain than the input image (due to the boundary condition
% adopted for the operator A), the output image must be cropped to be
% compared to u0 (indeed, the values computed outside of the domain of u
% are not relevant since they are only driven by the STV term, because no
% observation is availaible at those positions).  
[u,E] = solver_tv_relax(u0,lambda,A,adj_A,L_A,niter,'init',init,'verbose',verbose,'video',video,'displayFcn',crop_roi); 
u_crop = crop_roi(u); 
figure('Name','restored image'); imshow2(u_crop); 
figure(); plot(1:length(E),E); title('energy evolution'); xlabel('iteration'); ylabel('energy'); 

%% Play with the video mode
%
% set the optional parameter 'video' to 'true' to enable the video mode.
%
% when enabled, the image 'displayFcn(u)' is displayed at each iteration of
% the algorithm. 
%
% 'displayFcn' is also an optional input (Matlab function or Macro, with
% default value @(u)u).
%
% When video mode is enabled, the algorithms iterations can be stopped at
% any time by pressing key the key 'q' of the keyboard. 

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(crop_roi(u))))); % display log(1+abs(fftshift(fft2(crop_roi(u))))) at each iteration
solver_tv_relax(u0,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(crop_roi(u)),normalize(logdft(u))]; 
solver_tv_relax(u0,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',displayFcn); 
