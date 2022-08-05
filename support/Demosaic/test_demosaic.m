clearvars; close all;

addpath(fullfile('..','Validation'));
addpath(fullfile('..','Mosaic'));
addpath(fullfile('..','Load_Info'));
addpath(fullfile('..','Visualization'));

% I=double(imread('peppers.png'));
% I=I/max(I(:));
% DynamicRange=2^nextpow2(max(I(:))+1)-1;
% I=Load_Dataset_Pansharpening('Washington_cut256_4','request',{'GT'});
% I_GT=I{1};
I_GT=load_natural('CAVE1_RGB');

mask=load_mask('sizes',I_GT.size,'type','Bayer');
I_masked=mask.mosaic(I_GT,mask);

% methods={'RI','MLRI','MLRI2','ARI','ARI2'};
methods={'ID','IID','SD','ISD'};

I_out=demosaic_classic(I_masked,mask,methods);

[MR,qi_label,order]= indexes_evaluation_RR(I_out,I_GT,'qindex_list',{'ssim','psnr'});

figure;
subplot(221); a=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display)); imshow(a,[]); title('Original');
subplot(222); a=visualize_mask(mask,0); imshow(a,[]); title('Mask');
subplot(223); a=viewimage_outputonly(I_masked,[0,I_GT.DynamicRange]); imshow(a,[]); title('Mosaic');
subplot(224); a=viewimage_outputonly(I_out.data(:,:,I_GT.Bands_to_display,order(order(1,:)==1))); imshow(a,[]); title('Best Demosaic');