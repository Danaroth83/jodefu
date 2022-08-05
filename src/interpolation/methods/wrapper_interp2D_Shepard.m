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
% I_out=interp2D_Shepard(I_in,mask,Nsamples,p);
% 
% Input:
% I_In: Image with missing samples (sizes: Nx x Ny x Nb)
% mask: Position of missing samples marked with zeros (same sizes as I_in)
% Nsamples: Number of samples to be considered in the surroundings
%           (default: amount of non zero samples)
% p: Order of the interpolation
%
% Output:
% I_out: Reconstructed image (If the image is monocromatic it has 
%        sizes: Nx x Ny x Ns, otherwise Nx x Ny x Nb x Ns)
%

function I_out=wrapper_interp2D_Shepard(I_in,mask,Nsamples,p)

if nargin<=1 || isempty(mask), mask=ones(size(I_in)); mask(I_in==0)=0; end
if nargin<=2 || isempty(Nsamples), Nsamples=prod(size(I_in,1),size(I_in,2)); end
if nargin<=3 || isempty(p), p=2; end

%% Reshaping input

% Fix for 4D+ images
sizes=size(mask);
I_in=reshape(I_in,size(mask,1),size(mask,2),[]);
mask=reshape(mask,size(mask,1),size(mask,2),[]);
I_in=I_in.*mask; % Fix for mosaicked images

[y,x]=meshgrid(1:size(mask,2),1:size(mask,1));
x=x(:); y=y(:);

I_out=zeros([size(mask,1)*size(mask,2),size(mask,3),length(p)]);
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
    
    Nsamples_current=min(Nsamples,length(xi));
    fa=zeros(length(xf),length(p));
    
    for jj=1:length(xf)
        
        fprintf('sample: %d\n',jj);
        rf_ss=sqrt((xi-xf(jj)).^2+(yi-yf(jj)).^2);
        [rf_ss,idx_ss]=mink(rf_ss,Nsamples_current);
        fi_ss=fi(idx_ss);
        
        for ii=1:length(p)
            inv_rf_ss=1./rf_ss.^p(ii);
            fa(jj,ii)=sum(fi_ss.*inv_rf_ss)./sum(inv_rf_ss);
        end

    end
    I_out(:,kk,:)=repmat(I_current,[1,length(p)]);
    I_out(mask_zeroidx,kk,:)=fa;
end
I_out=reshape(I_out,[size(mask,1),size(mask,2),size(mask,3),length(p)]);
I_out=reshape(I_out,[sizes,length(p)]);
