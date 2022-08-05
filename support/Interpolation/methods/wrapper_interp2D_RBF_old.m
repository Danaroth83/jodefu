% RADIAL BASIS FUNCTION 2D SCATTERED INTERPOLATION
%
% Description:
% This script works as frontend for the toolbox of the Radial Basis
% Function (RBF) [1], applied for 2D interpolation, for images in a grid
% with missing samples, marked by the zeros in the mask. The script also
% supports multidimensional images, but the interpolation is just performed
% over each 
%
% Usage:
% I_out=RBF_interp2D(I_in,mask,Nsamples,sigma,type);
% Eg:
% im_orig=double(imread('peppers.png')); 
% im_orig=imresize(im_orig(:,:,1)./max(im_orig(:)),1/4);
% mask_gen=randn(size(im_orig)); [~,mask_idx]=max(mask_gen,[],3);
% mask=[]; for ii=1:size(mask_gen,3), mask=cat(3,mask,mask_idx==ii); end
% I_in=sum(im_orig.*mask,3);
% I_out=RBF_interp2D(I_in,mask,[],0.7,'Gaussian');
% subplot(221); imshow(im_orig); title('Original');
% subplot(222); imshow(mask); title('Mask');
% subplot(223); imshow(I_in); title('Mosaic Image');
% subplot(224); imshow(I_out); title('Reconstructed');
% 
% Input:
% I_In: Image with missing samples (sizes: Nx x Ny x Nb)
% mask: Position of missing samples (same sizes as I_in )
% Nsamples: Number of samples to be considered in the RBF surroundings
%           (default: amount of non zero samples)
% sigma: Parameter of the RBF; it's the reciprocal of the variance for a
%        Gaussian RBF (default: 1; can be inputted as vector of sizes 1 x Ns)
% type: The shape of the RBF (default: 'Gaussian')
%    -'Gaussian':              RBF(r)=exp(-sigma^2 r^2)
%    -'InverseQuadratic':      RBF(r)=1/(1+sigma^2 r^2)
%    -'MultiQuadratic':        RBF(r)=sqrt(1+sigma^2 r^2)
%    -'InverseMultiQuadratic': RBF(r)=1/sqrt(1+sigma^2 r^2)
%    -'ThinPlate':             RBF(r)=sigma^2*r^2*ln(sigma*r)
%
% Output:
% I_out: Reconstructed image (If the image is monocromatic it has 
%        sizes: Nx x Ny x Ns, otherwise Nx x Ny x Nb x Ns)
%
% Reference:
% [1] Sarra S., "The Matlab Radial Basis Function Toolbox". Journal of 
% Open Research Software, 5: 8,

function I_out=wrapper_interp2D_RBF_old(I_in,mask,~,sv,type)

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'RBF'));

if nargin<=1 || isempty(mask), mask=ones(size(I_in)); mask(I_in==0)=0; end
%if nargin<=2 || isempty(Nsamples), Nsamples=prod(size(image,1),size(image,2)); end
if nargin<=3 || isempty(sv), sv=1; end
if nargin<=4 || isempty(type), type='Gaussian'; end

%% Loading RBF

if any(strcmpi(type,{'Gaussian','gax'}))
    phi = gax();
elseif any(strcmpi(type,{'InverseQuadratic','iqx'}))
    phi = iqx();
elseif any(strcmpi(type,{'InverseMultiQuadratic','imqx'}))
    phi = imqx();
elseif any(strcmpi(type,{'MultiQuadratic','mqx'}))
    phi = mqx();
elseif any(strcmpi(type,{'ThinPlate','tpx'}))
    phi = tpx();
else
    error('Unknown RBF')
end

%% Parameters of the RBF solver
mu = 1.5e-13;
safe = false; % If false, it uses Cholesky factorization for regularization
% safe = true;

%% Reshaping input

% Fix for 4D+ images
sizes=size(mask);
I_in=reshape(I_in,size(mask,1),size(mask,2),[]);
mask=reshape(mask,size(mask,1),size(mask,2),[]);
I_in=I_in.*mask; % Fix for mosaicked images

[y,x]=meshgrid(1:size(mask,2),1:size(mask,1));
x=x(:); y=y(:);

I_out=zeros([size(mask,1)*size(mask,2),size(mask,3),length(sv)]);
tol=0.0001*max(mask(:));

for kk=1:size(I_in,3)

    %% Find positions of zeros
    I_current=reshape(I_in(:,:,kk),[],1);
    mask_current=reshape(mask(:,:,kk),[],1);
    mask_zeroidx=mask_current(:)<tol;

    %% Original samples
    xi=x(~mask_zeroidx);
    yi=y(~mask_zeroidx);
    fi=I_current(~mask_zeroidx)./(mask_current(~mask_zeroidx)).^2;

    %% Positions to interpolate
    xf=x(mask_zeroidx);
    yf=y(mask_zeroidx);
    
    %% Calculate distances between samples and target points
    ri = phi.distanceMatrix2d(xi,yi);
    rf = phi.distanceMatrix2d(xi,yi,xf,yf);

    init=zeros(size(I_in,1)*size(I_in,2),1);
    init(~mask_zeroidx)=fi;
    for ii=1:length(sv)
        s = sv(ii); % fprintf('s=%.1f\n',s);      
        B = phi.rbf(ri,s);
        a = phi.solve(B,fi,mu,safe);
        H = phi.rbf(rf,s);
        fa = H*a;
        I_recon=init;
        I_recon(mask_zeroidx)=fa;
        I_out(:,kk,ii)=I_recon;
        % er(ii) = sum((fa - f_ref).^2);
        % if ii==1 || er(ii)<er(ii-1), fa_save=fa; end
    end
end
I_out=reshape(I_out,[size(mask,1),size(mask,2),size(mask,3),length(sv)]);
I_out=reshape(I_out,[sizes,length(sv)]);



