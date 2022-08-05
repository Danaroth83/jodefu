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
%    -'Gaussian':              RBF(r)=exp(-(sigma*r)^2)
%    -'InverseQuadratic':      RBF(r)=1/(1+(sigma*r)^2)
%    -'MultiQuadratic':        RBF(r)=sqrt(1+(sigma*r)^2)
%    -'InverseMultiQuadratic': RBF(r)=1/sqrt(1+(sigma*r)^2)
%    -'ThinPlate':             RBF(r)=(sigma*r)^2*ln(sigma*r)
%
% Output:
% I_out: Reconstructed image (If the image is monocromatic it has 
%        sizes: Nx x Ny x Ns, otherwise Nx x Ny x Nb x Ns)
%
% Reference:
% [1] Sarra S., "The Matlab Radial Basis Function Toolbox". Journal of 
% Open Research Software, 5: 8,

function I_out=wrapper_interp2D_RBF_reject(I_in,mask,Nsamples,sv,type)

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'RBF'));
addpath(fullfile(current_folder,'knn'));

if nargin<=1 || isempty(mask), mask=ones(size(I_in)); mask(I_in==0)=0; end
if nargin<=2 || isempty(Nsamples), Nsamples=prod(size(I_in,1),size(I_in,2)); end
if nargin<=3 || isempty(sv), sv=1; end
if nargin<=4 || isempty(type), type='Gaussian'; end

%% Loading RBF

if any(strcmpi(type,{'Gaussian','gax'}))
    phi = gax();
    mu = 1.5e-13;
    safe = false;
elseif any(strcmpi(type,{'InverseQuadratic','iqx'}))
    phi = iqx();
    mu = 1.5e-13;
    safe = false;
elseif any(strcmpi(type,{'InverseMultiQuadratic','imqx'}))
    phi = imqx();
    mu = 1.5e-13;
    safe = false;
elseif any(strcmpi(type,{'MultiQuadratic','mqx'}))
    phi = mqx();
    mu = 1.5e-13;
    safe = false;
elseif any(strcmpi(type,{'ThinPlate','tpx'}))
    phi = tpx();
    mu = 1.5e-13;
    safe = true; % Matrix is not semidefinite positive
else
    error('Unknown RBF')
end


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
    % ri = phi.distanceMatrix2d(xi,yi);
    % rf = phi.distanceMatrix2d(xi,yi,xf,yf);
    
    Nsamples_current=min(Nsamples,length(xi));
    % rf=zeros(length(xf),Nsamples_current);
    % idx=zeros(length(xf),Nsamples_current);
    fa=zeros(length(xf),length(sv));
    
    % idx_ss=(nearestneighbour([xf,yf].',[xi,yi].', 'NumberOfNeighbours', Nsamples_current)).';
    % xi_ss=xi(idx_ss); yi_ss=yi(idx_ss); fi_ss=fi(idx_ss);
    % rf_ss=sqrt((xi_ss-xf).^2+(yi_ss-yf).^2);
    
    % idx_ss=knnsearch2([xf,yf],[xi,yi],Nsamples_current);
    % xi_ss=xi(idx_ss); yi_ss=yi(idx_ss); fi_ss=fi(idx_ss);
    % rf_ss=sqrt((xi_ss-xf).^2+(yi_ss-yf).^2);
    
    % rf_ss=pdist2([xf,yf],[xi,yi]);
    % [rf_ss,idx_ss]=mink(rf_ss,Nsamples_current,2);
    % xi_ss=xi(idx_ss); yi_ss=yi(idx_ss); fi_ss=fi(idx_ss);
    
    idx_ss=knnsearch([xf,yf],[xi,yi],Nsamples_current);
    xi_ss=xi(idx_ss); yi_ss=yi(idx_ss); fi_ss=fi(idx_ss);
    rf_ss=sqrt((xi_ss-xf).^2+(yi_ss-yf).^2);
    
    for jj=1:length(xf)
        ri_current = phi.distanceMatrix2d(xi_ss(jj,:),yi_ss(jj,:));
        fi_current = fi_ss(jj,:).';
        for ii=1:length(sv)
            s = sv(ii); % fprintf('s=%.1f\n',s);      
            B = phi.rbf(ri_current,s);
            a = phi.solve(B,fi_current,mu,safe);
            H = phi.rbf(rf_ss(jj,:),s);
            fa(jj,ii) = H*a;
        end
    end
    I_out(:,kk,:)=repmat(I_current,[1,length(sv)]);
    I_out(mask_zeroidx,kk,:)=fa;
end
I_out=reshape(I_out,[size(mask,1),size(mask,2),size(mask,3),length(sv)]);
I_out=reshape(I_out,[sizes,length(sv)]);



