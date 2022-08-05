%%         Discrete Total Variation based image zooming
% 
% This example illustrates how the generic 'solver_tv_constr.m' function
% can be used to perform image inpainting. 
%
% Given an input image u0, compute
%
%      argmin_{u in R^{nx*ny}} TV(u) subject to A(u) = u0
% 
% where A denotes a captor integration operator (which averages the values
% of the image over squared cells). 
%
% This script makes use of data (reference images) stored into the 'data'
% folder of this program, the following relative path must be adapted
% accordingly.  
%
relative_path_to_data = './data/'; 

%% Image zooming

ref = double(imread([relative_path_to_data,'eye.tif'])); % reference image 
z = 4; 
A = @(u)cell_averaging(u,z); 
u0 = A(ref); 
proj_constr = @(u) u + kron((u0-A(u)),ones(z,z));

figure('Name','reference image'); imshow2(ref); 
figure('Name','unzoomed'); imshow2(u0);

% set parameters, run algorithm and display result
[ny,nx] = size(ref); 
niter = 500;
verbose = true; 
video = true; 

[u,E] = solver_tv_constr(nx,ny,proj_constr,niter,'verbose',verbose,'video',video); 
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
solver_tv_constr(nx,ny,proj_constr,niter,'verbose',verbose,'video',video); 

% display only a particular region-of-interest
crop_roi = @(u)u(100:220,100:270); 
solver_tv_constr(nx,ny,proj_constr,niter,'verbose',verbose,'video',video,'displayFcn',crop_roi);

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
solver_tv_constr(nx,ny,proj_constr,niter,'verbose',verbose,'video',video,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
solver_tv_constr(nx,ny,proj_constr,niter,'verbose',verbose,'video',video,'displayFcn',displayFcn); 
