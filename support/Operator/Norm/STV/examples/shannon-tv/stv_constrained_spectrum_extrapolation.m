%%       Shannon Total Variation based spectrum extrapolation
% 
% This example illustrates how the generic 'solver_stv_constr.m' function
% can be used to perform image inpainting. 
%
% Given an input image u0, compute
%
%      argmin_{u in R^{nx*ny}} STVn(u) subject to A(u) = u0
% 
% where A denotes a frequency masking operator of the type 
%
% A(u) = ifft2(fft2(u).*mask), 
%
% mask being a symmetric binary mask given in the Fourier domain.   
%
% This script makes use of data (reference images) stored into the 'data'
% folder of this program, the following relative path must be adapted
% accordingly.  
%
relative_path_to_data = './data/'; 

%% Spectrum extrapolation 

% load & display a reference image and its spectrum amplitude (in log scale)
ref = double(imread([relative_path_to_data,'clock2.tif'])); ref = ref(288:end,1:225); % reference image 
figure('Name','reference image'); imshow2(ref); 
figure('Name','reference spectrum (log scale)'); imshow2(log(1+abs(fftshift(fft2(ref))))); 

% compute & display a rectangular frequency mask (or a radon frequency mask)
[ny,nx] = size(ref); 
mask = rectangularmask(nx,ny,round(nx/3),round(ny/3)); % rectangular mask
% mask = radonmask(nx,ny,64); % radon mask with 64 rays

% compute u0 = A(u)
A = @(u)real(ifft2(fft2(u).*mask)); 
u0 = A(ref); 
figure('Name','input spectrum'); imshow2(log(1+abs(fftshift(fft2(u0))))); 

% define the projection over the constraint set C = {u in R^{nx*ny} ,A(u)=u0}
dft_u0 = fft2(u0); 
proj_constr = @(u) real(ifft2(fft2(u).*(1-mask) + dft_u0.*mask));

% set parameters, run algorithm and display results 
[ny,nx] = size(ref); 
n = 2; 
niter = 300;
gain = 10*mean(mask(:)); 
verbose = true; 
video = true; 
[u,E] = solver_stv_constr(nx,ny,n,proj_constr,niter,'gain',gain,'verbose',verbose,'video',video); 

figure('Name','restored image'); imshow2(u); 
figure('Name','extrapolated spectrum'); imshow2(log(1+abs(fftshift(fft2(u))))); 
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
solver_stv_constr(nx,ny,n,proj_constr,niter,'gain',gain,'verbose',verbose,'video',video); 

% display only a particular region-of-interest
crop_roi = @(u)u(50:200,50:200);  
solver_stv_constr(nx,ny,n,proj_constr,niter,'gain',gain,'verbose',verbose,'video',video,'displayFcn',crop_roi); 

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
solver_stv_constr(nx,ny,n,proj_constr,niter,'gain',gain,'verbose',verbose,'video',video,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
solver_stv_constr(nx,ny,n,proj_constr,niter,'gain',gain,'verbose',verbose,'video',video,'displayFcn',displayFcn); 
