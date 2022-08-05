clearvars; close all;

current_folder=fileparts(mfilename('fullpath'));
support_folder=fullfile(current_folder,'..','..','support');
addpath(fullfile(support_folder,'Visualization'),...
        fullfile(support_folder,'Mosaic'));

input_folder=fullfile(current_folder,'..','..','output','test_compression');
output_folder=fullfile(current_folder,'..','..','output','figures_PhD','direct');

mkdir(output_folder);
im_tag='Beijing_cut256_WV3_WV3_4';
ratio=2;
ratio_scale=4;
ratio_mask=4;

request={'GT','EXP','MS_LR','PAN'};
I_in=Load_Dataset_Pansharpening(im_tag,'request',request,'ratio',ratio);
I_GT=I_in{1}; I_EXP=I_in{2}; I_MS_LR=I_in{3}; I_PAN=I_in{4};

for jj=2:3
    I_in{jj}.data=imresize(imresize(I_in{jj}.data,1/ratio_scale),ratio_scale);
end

color={[0,0,1],[0,1,0],[1,0,0],[1,1,0]};
color_label={'B','G','R','N'};
%y=opA.mosaic(cat(3,I_PAN.data,I_EXP.data));

%I_PAN_DEMO=opA.mosaic(cat(3,I_PAN.data,I_EXP.data));

label=request(1:3);
I_view=cell(1,numel(label));

for jj=1:numel(label)
    for kk=1:size(I_in{1}.data,3)
        a=viewimage_outputonly(I_in{jj}.data(:,:,kk));
        I_current=a(:,:,1);
        I_view{jj}=cat(3,I_view{jj},I_current);
        I_current=I_current.*shiftdim(color{kk},-1);
        %figure; imshow(I_current);
        %if kk==1 || kk==4, ratio_scale_curr=1; else, ratio_scale_curr=ratio_scale; end
        imwrite(I_current,fullfile(output_folder,[im_tag,'_',label{jj},'_',color_label{kk},'.png']));
    end
end
KerBlu=load_KerBlu( im_tag, 4, I_EXP.GNyq,min(size(I_EXP.data,1),size(I_EXP.data,2)),false);
for kk=1:size(I_MS_LR.data,3)
    a=KerBlu(1+15:end-15,1+15:end-15,kk);
    a=a/max(a(:));
    imwrite(imresize(ind2rgb(uint8(a*255),jet(255)),15,'nearest'),fullfile(output_folder,[im_tag,'_KERBLU_',color_label{kk},'.png']));
end

I_PAN_view=viewimage_outputonly(I_PAN.data);
imwrite(I_PAN_view,fullfile(output_folder,[im_tag,'_PAN.png']));
I_PAN_ds=imresize(I_PAN_view(:,:,1:size(I_PAN.data,3)),1/ratio_mask);
I_MS_LR_ds=imresize(I_view{3},1/ratio_mask);
I_EXP_ds=imresize(I_view{2},1/ratio_mask);

mask=load_mask('sizes',size(I_EXP_ds),'band_HR',size(I_in{4}.data,3),...
    'ratio',ratio,'mask_PAN','coverage','mask','period');
%opA =load_degmaskoperator('lpfilter',I_in{2}.KerBlu,'spectralweights',...
%    I_in{4}.spectralweights,'mask',mask,'flag_fft',0,'radius',1);
y=mask.mosaic(cat(3,I_PAN_ds,I_EXP_ds),mask);
x_demo=mask.demosaic(y,mask);
I_MS_LR_MOS=x_demo.image_LR;
I_mosaic_view=zeros([size(I_MS_LR_ds,1),size(I_MS_LR_MOS,2),3]);
I_msmask_view=zeros([size(I_MS_LR_ds,1),size(I_MS_LR_MOS,2),3]);
I_combo_view=repmat(y.*(mask.data(:,:,1)),[1,1,3]);

I_PAN_MOS=I_combo_view(:,:,1); I_PAN_MOS(mask.data(:,:,1)==0)=1;

imwrite(imresize(I_PAN_MOS,ratio_mask,'nearest'),fullfile(output_folder,[im_tag,'_PAN_MOSAIC.png']));
imwrite(imresize(1-x_demo.mask_HR,ratio_mask,'nearest'),fullfile(output_folder,[im_tag,'_PAN_MASK.png']));

for kk=1:size(I_in{1}.data,3)
    I_msmask_view=I_msmask_view+(x_demo.mask_LR(:,:,kk)).*shiftdim(color{kk},-1);
    I_mosaic_view=I_mosaic_view+I_MS_LR_MOS(:,:,kk).*shiftdim(color{kk},-1);
    I_combo_view=I_combo_view+y.*(mask.data(:,:,kk+1)).*shiftdim(color{kk},-1);
end
imwrite(imresize(I_mosaic_view,ratio_mask,'nearest'),fullfile(output_folder,[im_tag,'_MS_LR_MOSAIC.png']));
imwrite(imresize(I_msmask_view,ratio_mask,'nearest'),fullfile(output_folder,[im_tag,'_MS_LR_MASK.png']));
imwrite(imresize(I_combo_view,ratio_mask,'nearest'),fullfile(output_folder,[im_tag,'_MOSAIC.png']));

ra=1.5;
size_filter=21;  % length of the filter    
stopband=1.3; attenuation=30; % Parameters to design Butterworth filter
[f1,f2]=freqspace(201,'meshgrid');   
[n,Wn]=buttord(1/ra,min(stopband/ra,1-eps),10*log10(2),attenuation);
butt_filt=1./sqrt(1+((f1.^2+f2.^2)/Wn^2).^n);
h=fwind1(butt_filt,hamming(size_filter));
% h=fwind1(butt_filt,boxcar(size_filter));
h=h./sum(h(:));
a=(abs(h)./max(abs(h(:))));
imwrite(ind2rgb(uint8(imresize(a(1+6:end-6,1+6:end-6)*255,15,'nearest')),jet(100)),fullfile(output_folder,[im_tag,'_FILTER.png']));

ra=3;
size_filter=21;  % length of the filter    
stopband=1.3; attenuation=30; % Parameters to design Butterworth filter
[f1,f2]=freqspace(201,'meshgrid');   
[n,Wn]=buttord(1/ra,min(stopband/ra,1-eps),10*log10(2),attenuation);
butt_filt=1./sqrt(1+((f1.^2+f2.^2)/Wn^2).^n);
h=fwind1(butt_filt,hamming(size_filter));
% h=fwind1(butt_filt,boxcar(size_filter));
h=h./sum(h(:));

I_PAN_ra=imfilter(I_PAN_view,h);
imwrite(I_PAN_ra,fullfile(output_folder,[im_tag,'_PAN_ra.png']));