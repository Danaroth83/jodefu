% SHEPARD INTERPOLATION
%
% Description:
% This script implements the Inverse Weighting Interpolation (also known 
% as Shepard interpolation) for 2D images with missing samples, whose
% positions are marked in a mask. Multidimensional images are supported.
% Specifically Shepard interpolation for a sample y is defined as:
%    y=(sum_i w_i x_i) / (sum_i w_i), where: 
%     x_i is a given input sample
%     w_i = 1/r_i^p 
%     r_i is the distance between the missing sample and the position of x_i
%     p is the order of the interpolation
%
% Usage:
% I_out=interp2D_Shepard(xi,yi,fi,xo,yo,Nsamples,p,flag_modified);
% 
% Input:
% xi,yi: Column vectors for the input sample coordinates
% fi: Column vector for the input values
% xo,yo: Column vectors for the output coordinates
% Nsamples: Number of samples to be considered in the surroundings
%           (default: amount of non zero samples)
% p: Order of the interpolation (default: 2)
% method: if 'classic', calculates the weights as w_i=(1/r_i)^p
%         if 'mod', as w_i=(1/r_i-1/R)^p, where R is the max radius of
%                distance between a given output sample and the input samples
%         (default: 'mod')
%
% Output:
% I_out: Reconstructed image (If the image is monocromatic it has 
%        sizes: Nx x Ny x Ns, otherwise Nx x Ny x Nb x Ns)

function fo=interp2D_Shepard(xi,yi,fi,xo,yo,Nsamples,p,method)

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'MRBFT-1.0'));

if nargin<=5 || isempty(Nsamples), Nsamples=prod(size(I_in,1),size(I_in,2)); end
if nargin<=6 || isempty(p), p=2; end
if nargin<=7 || isempty(method), method='mod'; end

Nsamples=min(Nsamples,length(xi));
fo=zeros(length(xo),length(p));
fprintf('Sample:           ')
    
for jj=1:length(xo)

    fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d',jj,length(xo));
    rf_ss=sqrt((xi-xo(jj)).^2+(yi-yo(jj)).^2);
    [rf_ss,idx_ss]=mink(rf_ss,Nsamples);
    fi_ss=fi(idx_ss);

    if strncmpi(method,'mod',3), max_rf_ss=max(rf_ss); end
    for ii=1:length(p)
        if strncmpi(method,'mod',3)
            w=(1./rf_ss-1/max_rf_ss).^p(ii);
        else
            w=1./rf_ss.^p(ii);
        end
        fo(jj,ii)=sum(fi_ss.*w)./sum(w);
    end

end
fprintf('\n')
