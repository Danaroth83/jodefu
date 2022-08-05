%%       Shannon Total Variation based relaxed image inpainting
% 
% This example illustrates how the generic 'solver_stv_relax.m' function
% can be used to perform relaxed image inpainting. 
%
% Given an input image u0, compute a minimizer of the energy
%
%           E(u) := ||A(u)-u0||^2 +lambda*STVn(u)
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

%% Relaxed image inpainting

% load & display the reference image
ref = double(imread([relative_path_to_data,'fiat.tif'])); % reference image 
figure('Name','reference image'); imshow2(ref); 

% define operator A, its adjoint adj_A, and an upper bound L_A for |||A|||
mask = double(rand(size(ref)) > 0.5); 
A = @(u)u.*mask; % masking operator
adj_A = A; % A is self-adjoint ...
L_A = 1; % ... with unitary l2 induced norm (as soon as mask is non identically zero)

% compute & display the masked and noisy image u0 
u0 = A(ref) + 2*randn(size(ref)); 
figure('Name','masked and noisy image'); imshow2(u0);

% set parameters, run algorithm and display result
lambda = 0.5; 
n = 2; 
niter = 200;
init = adj_A(u0); 
gain = 10*mean(mask(:)); 
verbose = true; 
video = true; 
[u,E] = solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'gain',gain,'verbose',verbose,'video',video); 
figure('Name','restored image'); imshow2(u); 
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
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'gain',gain,'verbose',true,'video',true); 

% display only a particular region-of-interest
crop_roi = @(u)u(100+(1:200),100+(1:200)); 
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'gain',gain,'verbose',true,'video',true,'displayFcn',crop_roi); 

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'gain',gain,'verbose',true,'video',true,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'gain',gain,'verbose',true,'video',true,'displayFcn',displayFcn); 
