%% TEST 1: DETERMINISTIC MASK ASSIGNED TO A SINGLE BAND
clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig(:,:,1)./max(im_orig(:)),1/4);
mask=ones(size(im_orig));
mask(2:2:end,2:2:end)=0;
I_in=sum(im_orig.*mask,3);

% method='RBF';
method='RBF_spline';
% method='RBF_iq';
% method='RBF_mq';
% method='RBF_imq';
% method='WB';
% method='NN';
% method='IDW';

I_out=interp2D_mosaic(I_in,'mask',mask,'method',method,'Nsamples',100);

figure;
subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out); title('Reconstructed');

%% TEST 2: BAYER MASK

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask=zeros(size(im_orig));
mask(1:2:end,2:2:end,1)=1; mask(2:2:end,1:2:end,3)=1;
mask(1:2:end,1:2:end,2)=1; mask(2:2:end,2:2:end,2)=1;
I_in=sum(im_orig.*mask,3);

method='RBF_spline';
I_out=interp2D_mosaic(I_in,'mask',mask,'method',method,'Nsamples',100);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
subplot(224); imshow(I_out); title('Reconstructed');

%% TEST 3: RANDOM MASK

clearvars; close all;
im_orig=double(imread('peppers.png')); 
im_orig=imresize(im_orig./max(im_orig(:)),1/4);
mask_gen=rand(size(im_orig)); [~,mask_idx]=max(mask_gen,[],3);
mask=zeros(size(mask_gen));
for ii=1:3, mask(:,:,ii)=mask_idx==ii; end
I_in=sum(im_orig.*mask,3);

method='RBF_spline';
I_out=interp2D_mosaic(I_in,'mask',mask,'method',method,'Nsamples',100);

% MSE_ev=zeros(size(I_out,4),1);
% for ii=1:size(I_out,4), MSE_ev(ii)=sum((im_orig(:)-reshape(I_out(:,:,:,ii),[],1)).^2)/numel(I_in); end
% [~,idx_MSE]=min(MSE_ev);

subplot(221); imshow(im_orig); title('Original');
subplot(222); imshow(mask); title('Mask');
subplot(223); imshow(I_in); title('Mosaic Image');
% subplot(224); imshow(I_out(:,:,:,idx_MSE)); title('Reconstructed');
subplot(224); imshow(I_out); title('Reconstructed');

