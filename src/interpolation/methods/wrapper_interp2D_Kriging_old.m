% KRIGING 2D SCATTERED INTERPOLATION
%
% Description:
% This script works as frontend for the toolbox of the Radial Basis
% Function (RBF) [1], applied for 2D interpolation, for images in a grid
% with missing samples, marked by the zeros in the mask. The script also
% supports multidimensional images, but the interpolation is just performed
% over each 
%
% Usage:
% I_out=RBF_Kriging(I_in,mask,Nsamples,Nbins,maxdist);
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

function I_out=interp2D_Kriging_old(I_in,mask,Nsamples,model,Nbins)

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'MRBFT-1.0'));
addpath(fullfile(current_folder,'Kriging'));

if nargin<=1 || isempty(mask), mask=ones(size(I_in)); mask(I_in==0)=0; end
if nargin<=2 || isempty(Nsamples), Nsamples=prod(size(I_in,1),size(I_in,2)); end
if nargin<=3 || isempty(model), model=[]; end
if nargin<=4 || isempty(Nbins), Nbins=20; end


%% Reshaping input

% Fix for 4D+ images
sizes=size(mask);
I_in=reshape(I_in,size(mask,1),size(mask,2),[]);
mask=reshape(mask,size(mask,1),size(mask,2),[]);
I_in=I_in.*mask; % Fix for mosaicked images

[y,x]=meshgrid(1:size(mask,2),1:size(mask,1));
x=x(:); y=y(:);

I_out=zeros([size(mask,1)*size(mask,2),size(mask,3)]);
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
    
    %% Calculate variogram
    
    Nsamples_current=min(Nsamples,length(xi));
    fa=zeros(length(xf),1);
    rfsq_ss=(xi.'-xf).^2+(yi.'-yf).^2;
    idx_ss=mink(rfsq_ss,Nsamples_current,2);
    xi_ss=xi(idx_ss); yi_ss=yi(idx_ss); fi_ss=fi(idx_ss); 
    
    % vario=variogram([xi,yi],fi,'nrbins',Nbins);
    
    ri_corr=reshape(sqrt((xi-xi.').^2+(yi-yi.').^2),[],1);
    limit=[];
    if isempty(limit), limit=max(ri_corr(:))/2; end
    % if isempty(limit), limit=min(max(xi(:))-min(xi(:)),max(yi(:))-min(yi(:)))/2; end
    idx_limit = (ri_corr<=limit);
    fi_corr=reshape(fi-fi.',[],1);
    ri_lim=ri_corr(idx_limit);
    fi_lim=fi_corr(idx_limit);
    [n,edges,bin]=histcounts(ri_lim,Nbins); % edges= limits of the bins, n= amount of elements in the bin, bin=index for elements of the bin
    value=zeros([Nbins,1]);
    for ii=1:Nbins
        value(ii)=sum(fi_lim(bin==ii).^2)/n(ii)/2;
    end
   vario.val=value; vario.distance=(edges(2:end)++edges(1:end-1)).'/2; vario.count=n.';
    
    [~,~,~,struct_vario]=variogramfit(vario.distance,vario.val,[],[],[],'model',model);
    
    fprintf('Current sample:             ');
    for jj=1:length(xf)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d',jj,length(xf));
        % fa(jj)=kriging(struct_vario,xi_ss(jj,:).',yi_ss(jj,:).',fi_ss(jj,:).',xf(jj),yf(jj));
        xi_current=xi_ss(jj,:).'; yi_current=yi_ss(jj,:).'; fi_current=fi_ss(jj,:).';
        % var_gamma=(fi_current-fi_current.').^2./sqrt((xi_current-xi_current.').^2+(yi_current-yi_current.').^2)/2;
        % var_gamma=cov(fi_current);
        var_gamma=struct_vario.func([struct_vario.range,struct_vario.sill],sqrt((xi_current-xi_current.').^2+(xi_current-xi_current.').^2));
        if ~isempty(struct_vario.nugget), var_gamma=var_gamma+struct_vario.nugget; end
        var_gamma=padarray(var_gamma,[1,1],1,'post'); var_gamma(end,end)=0;
        cov_gamma=struct_vario.func([struct_vario.range,struct_vario.sill],sqrt((xi_current-xf(jj)).^2+(yi_current-yf(jj)).^2));
        if ~isempty(struct_vario.nugget), cov_gamma=cov_gamma+struct_vario.nugget; end
        cov_gamma=cat(1,cov_gamma,1);
        w=pinv(var_gamma)*cov_gamma; w=w(1:end-1);
        fa(jj)=sum(w.*fi_current);
    end
    fprintf('. Done!\n')
    
    
%     for jj=1:length(xf)
%         ri_current = phi.distanceMatrix2d(xi_ss(jj,:),yi_ss(jj,:));
%         fi_current = fi_ss(jj,:).';
%         for ii=1:length(sv)
%             s = sv(ii); % fprintf('s=%.1f\n',s);      
%             B = phi.rbf(ri_current,s);
%             a = phi.solve(B,fi_current,mu,safe);
%             H = phi.rbf(rf_ss(jj,:),s);
%             fa(jj,ii) = H*a;
%         end
%     end
    I_out(:,kk)=I_current;
    I_out(mask_zeroidx,kk)=fa;
end
I_out=reshape(I_out,[size(mask,1),size(mask,2),size(mask,3)]);
I_out=reshape(I_out,sizes);



