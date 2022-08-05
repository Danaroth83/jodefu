%%       Shannon Total Variation based relaxed image zooming
% 
% This example illustrates how the generic 'solver_stv_relax.m' function
% can be used to perform relaxed spectrum extrapolation. 
%
% Given an input image u0, compute a minimizer of the energy
%
%           E(u) := ||A(u)-u0||^2 +lambda*STVn(u)
% 
% where A denotes a captor integration operator (which averages the values
% of the image over squared cells). 
%
% This script makes use of data (reference images) stored into the 'data'
% folder of this program, the following relative path must be adapted
% accordingly.  
%
relative_path_to_data = './data/'; 

%% Relaxed image zooming

% define operator A, its adjoint adj_A, and an upper bound L_A for |||A|||
z = 4; % side lengths of the (square) integration cell 
A = @(u)cell_averaging(u,z); % cell averaging operator
adj_A = @(v)kron(v,ones(z,z)/(z^2)); % adjoint of A
L_A = 1/z; % upper bound for the l2 induced norm of A

% load & display the reference image and its unzoomed & noisy version
ref = double(imread([relative_path_to_data,'eye.tif'])); % reference image 
u0 = A(ref)+2*randn(size(A(ref))); % unzoomed and noisy image
figure('Name','reference image'); imshow2(ref); 
figure('Name','unzoomed and noisy'); imshow2(u0);

% set parameters, run algorithm and display result
lambda = 0.05; 
n = 2; 
niter = 100;
init = adj_A(u0)*z^2; 
verbose = true; 
video = true; 
[u,E] = solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',verbose,'video',video); 
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
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true); 

% display only a particular region-of-interest
crop_roi = @(u)u(100:220,100:270); 
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',crop_roi); 

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
solver_stv_relax(u0,n,lambda,A,adj_A,L_A,niter,'init',init,'verbose',true,'video',true,'displayFcn',displayFcn); 

