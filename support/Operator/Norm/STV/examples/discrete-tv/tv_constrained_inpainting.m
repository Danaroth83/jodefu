%%       Discrete Total Variation based image inpainting
% 
% This example illustrates how the generic 'solver_tv_constr.m' function
% can be used to perform image inpainting. 
%
% Given an input image u0, compute
%
%      argmin_{u in R^{nx*ny}} TV(u) subject to A(u) = u0
% 
% where A denotes a masking operator of the type 
%
% A(u) = u.*mask,
%
% mask being a binary mask.   
%
% This script makes use of data (reference images) stored into the 'data'
% folder of this program, the following relative path must be adapted
% accordingly.  
%
relative_path_to_data = './data/'; 

%% Image inpainting

% load a reference image, compute a random mask, and compute u0 = A(ref)
ref = double(imread([relative_path_to_data,'fiat.tif'])); % reference image 
mask = double(rand(size(ref)) > 0.6); % aroud 60 % of the mask is 0
A = @(u) u.*mask; 
u0 = A(ref); 

figure('Name','reference image'); imshow2(ref); 
figure('Name','binary mask'); imshow2(mask); 
figure('Name','input image'); imshow2(u0); 

% define the projection over the constraint set C = {u in R^{nx*ny} ,A(u)=u0}
proj_constr = @(u) u.*(1-mask) + u0.*mask; 

% set parameters, run algorithm and display results 
[ny,nx] = size(ref); 
niter = 300;
gain = 100*mean(mask(:)); 
verbose = true; 
video = true; 
figure('Name','reference image'); imshow2(ref); 
figure('Name','binary mask'); imshow2(mask); 
figure('Name','input image'); imshow2(u0); 
[u,E] = solver_tv_constr(nx,ny,proj_constr,niter,'gain',gain,'verbose',verbose,'video',video); 

figure('Name','restored image'); imshow2(u); 
figure(); plot(1:length(E),E); title('energy evolution'); xlabel('iteration'); ylabel('energy'); 

% check that ouptut image satisfies the constraint 
fprintf('computed ||A(u)-u0||_inf = %.10g\n',max(max(abs(A(u)-u0)))); 

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

% use video mode with default 'displayFcn': display u at each iteration
solver_tv_constr(nx,ny,proj_constr,niter,'gain',gain,'verbose',true,'video',true); 

% display only a particular region-of-interest
crop_roi = @(u)u(100+(1:200),100+(1:200)); 
solver_tv_constr(nx,ny,proj_constr,niter,'gain',gain,'verbose',true,'video',true,'displayFcn',crop_roi); 

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
solver_tv_constr(nx,ny,proj_constr,niter,'gain',gain,'verbose',true,'video',true,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
solver_tv_constr(nx,ny,proj_constr,niter,'gain',gain,'verbose',true,'video',true,'displayFcn',displayFcn); 

