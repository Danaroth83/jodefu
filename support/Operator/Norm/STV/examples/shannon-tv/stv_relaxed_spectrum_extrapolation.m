%%     Shannon Total Variation based relaxed spectrum extrapolation
% 
% This example illustrates how the generic 'solver_stv_relax.m' function
% can be used to perform relaxed spectrum extrapolation. 
%
% Given an input image u0, compute a minimizer of the energy
%
%           E(u) := ||A(u)-u0||^2 +lambda*STVn(u)
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

%% Relaxed spectrum extrapolation 

% load & display a reference image and its spectrum amplitude (in log scale)
ref = double(imread([relative_path_to_data,'clock2.tif'])); ref = ref(288:end,1:225); % reference image 
figure('Name','reference image'); imshow2(ref); 
figure('Name','reference spectrum (log scale)'); imshow2(log(1+abs(fftshift(fft2(ref))))); 

% compute & display a rectangular frequency mask (or a radon frequency mask)
[ny,nx] = size(ref); 
mask = rectangularmask(nx,ny,round(nx/3),round(ny/3)); % rectangular mask
% mask = radonmask(nx,ny,64); % radon mask with 64 rays
figure('Name','frequency mask'); imshow2(fftshift(mask)); 

% define operator A, its adjoint adj_A, and an upper bound L_A for |||A|||
A = @(u)real(ifft2(fft2(u).*mask)); 
adj_A = A; % A is self adjoint ...
L_A = 1; % ... with unitary l2 induced norm (as soon as mask is non identically zero)

% compute & display the input image u0 
u0 = A(ref)+0*randn(size(ref)); 
figure('Name','input image'); imshow2(u0);
figure('Name','input spectrum (log scale)'); imshow2(log(1+abs(fftshift(fft2(u0))))); 

% set parameters, run algorithm and display result
lambda = 0.4; 
n = 2; 
niter = 200;
gain = 10*mean(mask(:)); 
init = adj_A(u0); 
verbose = true; 
video = true; 
[u,E] = solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'gain',gain,'verbose',verbose,'video',video); 
figure('Name','restored image'); imshow2(u); 
figure('Name','extrapolated spectrum (log scale)'); imshow2(log(1+abs(fftshift(fft2(u))))); 
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

% use video mode with default 'displayFcn': display u at each iteration
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true); 

% display only a particular region-of-interest
crop_roi = @(u)u(50:200,50:200); 
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',crop_roi); 

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',displayFcn); 
