%%          Discrete Total Variation based image denoising
% 
% This example illustrates how the 'tvdenoise.m' function can be used to
% perform image denoising.  
%
% Given an input image u0, compute a minimizer of the energy
%
%           E(u) := ||u-u0||^2 +lambda*TVn(u)
%         
% which correspond to the classical Rudin-Osher-Fatemi model. 
%
% This script makes use of data (images) stored into the 'data' folder of
% this program, the following relative path must be adapted accordingly. 
%
relative_path_to_data = './data/'; 

%% Image denoising

% load & display the reference image and its blurry and noisy version
ref = double(imread([relative_path_to_data,'clock2.tif'])); % reference image
u0 = ref+10*randn(size(ref)); % blurry and noisy image
figure('Name','reference image'); imshow2(ref); 
figure('Name','blurry and noisy'); imshow2(u0); 

% set parameters, run algorithm, and display results
lambda = 7; 
niter = 100;
verbose = true; 
video = true; 
[u,E] = tvdenoise(u0,lambda,niter,'verbose',verbose,'video',video); 
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

% use video mode to display the spectrum evolution
logdft = @(u)log(1+abs(fftshift(fft2(u)))); % display log(1+abs(fftshift(fft2(u)))) at each iteration
tvdenoise(u0,lambda,niter,'verbose',true,'video',true,'displayFcn',logdft); 

% display both u and its spectrum
normalize = @(u) (u-min(u(:)))/(max(u(:))-min(u(:))); % contrast change
displayFcn = @(u) [normalize(u),normalize(logdft(u))]; 
tvdenoise(u0,lambda,niter,'verbose',true,'video',true,'displayFcn',displayFcn); 
