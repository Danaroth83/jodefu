%% TEST 1: DETERMINISTIC MASK ASSIGNED TO A SINGLE BAND
clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig(:,:,1)./max(im_orig(:)),1/4);
mask=ones(size(im_orig));
mask(2:2:end,2:2:end)=0;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF_old(I_in,mask,[],0.7,'Gaussian');
subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out); title('Reconstructed');

%% TEST 2: DETERMINISTIC MASK - 2 BANDS

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig(:,:,[1,2])./max(im_orig(:)),1/4);
mask=ones(size(im_orig));
mask(2:2:end,2:2:end,1)=0; mask(1:2:end,1:2:end,1)=0;
mask(2:2:end,1:2:end,2)=0; mask(1:2:end,2:2:end,2)=0;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF_old(I_in,mask,[],0.7,'Gaussian');
subplot(221); imshow(im_orig(:,:,[1,2,1])); title('Original');
subplot(222); imshow(mask(:,:,[1,2,1])); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,[1,2,1])); title('Reconstructed');

%% TEST 3: RANDOM MASK - 1 BAND

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig(:,:,1)./max(im_orig(:)),1/4);
mask=rand(size(im_orig));
mask(mask>1/4)=1; mask(mask<=1/4)=0;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF_old(I_in,mask,[],0.6,'Gaussian');
MSE_ev=zeros(size(I_out,3),1);
for ii=1:size(I_out,3), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);
subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,idx_MSE)); title('Reconstructed');

%% TEST 4: RANDOM MASK

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask_gen=rand(size(im_orig)); [~,mask_idx]=max(mask_gen,[],3);
mask=zeros(size(mask_gen));
for ii=1:3, mask(:,:,ii)=mask_idx==ii; end
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF_old(I_in,mask,[],0.5,'Gaussian');
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');

%% TEST 5: Bayer Mask

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask=zeros(size(im_orig));
mask(1:2:end,2:2:end,1)=1; mask(2:2:end,1:2:end,3)=1;
mask(1:2:end,1:2:end,2)=1; mask(2:2:end,2:2:end,2)=1;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF_old(I_in,mask,[],0.3:0.1:0.7,'Gaussian');
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');


%% TEST 6: LIMITED SAMPLES
clearvars; close all;
im_orig=double(imread('peppers.png'));
ratio=2;
im_orig=imresize(im_orig(:,:,1)./max(im_orig(:)),1/ratio);
mask=ones(size(im_orig));
mask(2:2:end,2:2:end)=0;
I_in=sum(im_orig.*mask,3);
% tic; I_out=interp2D_RBF_reject(I_in,mask,100,0.7,'Gaussian'); toc; 
tic; I_out=interp2D_RBF(I_in,mask,100,0.7,'Gaussian'); toc; 
subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out); title('Reconstructed');


%% TEST 7: LIMITED SAMPLES WITH BAYER

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask=zeros(size(im_orig));
mask(1:2:end,2:2:end,1)=1; mask(2:2:end,1:2:end,3)=1;
mask(1:2:end,1:2:end,2)=1; mask(2:2:end,2:2:end,2)=1;
I_in=sum(im_orig.*mask,3);
% tic; I_out=interp2D_RBF_old(I_in,mask,100,0.3:0.1:0.7,'Gaussian'); toc;
tic; I_out=interp2D_RBF(I_in,mask,100,0.3:0.1:0.7,'Gaussian'); toc;
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');

%% TEST 8: LIMITED SAMPLES FULL IMAGE (BAYER) VERY SLOW

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=im_orig./max(im_orig(:));
mask=zeros(size(im_orig));
mask(1:2:end,2:2:end,1)=1; mask(2:2:end,1:2:end,3)=1;
mask(1:2:end,1:2:end,2)=1; mask(2:2:end,2:2:end,2)=1;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF(I_in,mask,100,0.6,'Gaussian');
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');

%% TEST 9: SHEPARD LIMITED SAMPLES IMAGE (BAYER)

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask=zeros(size(im_orig));
mask(1:2:end,2:2:end,1)=1; mask(2:2:end,1:2:end,3)=1;
mask(1:2:end,1:2:end,2)=1; mask(2:2:end,2:2:end,2)=1;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_Shepard(I_in,mask,100,5:10);
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');

%% TEST 10: THIN PLATE SPLINE

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask=zeros(size(im_orig));
mask(1:2:end,2:2:end,1)=1; mask(2:2:end,1:2:end,3)=1;
mask(1:2:end,1:2:end,2)=1; mask(2:2:end,2:2:end,2)=1;
I_in=sum(im_orig.*mask,3);
% tic; I_out=interp2D_RBF_old(I_in,mask,100,0.3:0.1:0.7,'Gaussian'); toc;
tic; I_out=interp2D_RBF(I_in,mask,100,1,'ThinPlate'); toc;
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');

%% TEST 11: THIN PLATE SPLINE (RANDOM MASK)

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask_gen=rand(size(im_orig)); [~,mask_idx]=max(mask_gen,[],3);
mask=zeros(size(mask_gen));
for ii=1:3, mask(:,:,ii)=mask_idx==ii; end
I_in=sum(im_orig.*mask,3);
I_out=interp2D_RBF(I_in,mask,100,0.6:0.2:1.4,'ThinPlate');
MSE_ev=zeros(size(I_out,4),1);
for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
[~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');

%% TEST 12: DETERMINISTIC MASK ASSIGNED TO A SINGLE BAND (KRIGING)
clearvars; close all;
Nbins=40; model='spherical';
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig(:,:,1)./max(im_orig(:)),1/4);
mask=ones(size(im_orig));
mask(2:2:end,2:2:end)=0;
I_in=sum(im_orig.*mask,3);
I_out=interp2D_Kriging_old(I_in,mask,100,model,Nbins);
figure;
subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out); title('Reconstructed');

