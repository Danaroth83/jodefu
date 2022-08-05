function Enhanced_HSI = GFPCA(HSI, G_im, No_PCs, r, eps)
%GFPCA perform image fusion of hyperspectral and pan/RGB image, by using
%guided filter in PCA domain
%
%   [Enhanced_HSI,CTime] = GFPCA(HSI, G_im, No_PCs, sigma, eps)
%
%         Input:
%                 HSI                     - low-resolution hyperspectral image
%                 G_im                    - the guidance image, pan/RGB with high-resolution
%                 No_PCs                  - number of principal components used by guided filtering with pan/RGB image
%                 r                       - the size of local sliding window
%                 eps                     - regularization parameter determining the degree of blurring for the guided filter
%
%         Output:
%                 Enhanced_HSI            - enhanced hyperspectral image with both high spectral and spatial resolution
%
%
%       Copyright notes
%       Author: Wenzhi Liao IPI, Telin, Ghent University, Belgium
%       Email: wenzhi.liao@telin.ugent.be, wenzi.liao@gmail.com
%       Date: 30/10/2014
%
%    For those who want to use our codes, please cited our paper for the outcome of 2014 IEEE Data Fusion Contest, we are also glad to be co-authors.


G_im=double(G_im);
HSI=double(HSI);
[width,height, IRGB]=size(G_im);
[rows,cols,bands]=size(HSI);

%% decorrelate the HS image bands separating the information content from noise by PCA content from noise.
Xsp=reshape(HSI,rows*cols,bands);
W=pca2(Xsp,bands);
PCs=Xsp*W;
PCs=reshape(PCs,rows,cols,bands);
clear Xsp
  
%% Enhanced the first few principal components by pan/RGB image using guided filter
PCs_k= imresize(PCs(:,:,1:No_PCs),[width height], 'bicubic');
PCs_GF=zeros(width, height,No_PCs);
for i=1:No_PCs
    PCi=PCs_k(:,:,i);
    if IRGB==1
        PCs_GF(:,:,i)=guidedfilter(G_im,PCi,r,eps); % Fusion of HSI and pan image
    else
        PCs_GF(:,:,i)=guidedfilter_color(G_im,PCi,r,eps); % Fusion of HSI and RGB image
    end
end

%% Simply remove noise in the remaining principal components
PCs_B_k=zeros(width,height,bands-No_PCs);
for i=1:bands-No_PCs
      [thr,sorh,keepapp] = ddencmp('den','wv',PCs(:,:,i+No_PCs));
      PCs(:,:,i+No_PCs) = wdencmp('gbl',PCs(:,:,i+No_PCs),'sym4',4,thr,sorh,keepapp);
      PCs_B_k(:,:,i)=imresize(PCs(:,:,i+No_PCs),[width height], 'bicubic');
end


%% back project into original space
cleanPCA=cat(3,PCs_GF,PCs_B_k);
YPCA=reshape(cleanPCA,width*height,bands);YPCA=YPCA*W';
clear cleanPCA PCs_B_k PCs_GF
Enhanced_HSI=reshape(YPCA,width,height,bands);
clear YPCA

 
