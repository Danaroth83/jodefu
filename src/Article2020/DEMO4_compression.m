clearvars; close all;

support_folder=fullfile('..','..','support');
output_folder=fullfile('..','..','output','test_compression');
addpath(fullfile(support_folder,'Load_info'),...
        fullfile(support_folder,'Fusion'),...
        fullfile(support_folder,'Validation'),...
        fullfile(support_folder,'Visualization'));

% im_tag='Hobart1_RGB';
% im_tag='Hobart2_RGB';
% im_tag='Rome_RGB';
% im_tag='Sofia_cut2_ALI_ALI_RGB';
% im_tag='Hobart2';
% im_tag='China';
% im_tag='RdJ_WV3_WV3_4';
% im_tag='RdJ_WV3_WV3';
% im_tag='Rio_WV2_WV2';
% im_tag='Rome_WV2_WV2';
% im_tag='Washington_cut256_RGB';
im_tag='Washington_cut256_4';
% im_tag='Stockholm_cut256_RGB';
% im_tag='Stockholm_cut256_all';
% im_tag='Janeiro_cut256_all';
% im_tag='Janeiro_cut256_RGB';
% im_tag='SanFrancisco_QB_QB_4';
% im_tag='Beijing_WV3_WV3_4';
% im_tag='Hobart';

ratio=2;
flag_PANfromGT=0;
methods_list={'EXP','GSA','MTF-GLP-CBD','BDSD','BayesNaive'};
qindex_list={'SSIM','PSNR','ERGAS','SAM','SCC','UIQI','Q2^n'};


snr_PANfromGT=25;
I_cell=Load_Dataset_Pansharpening(im_tag,'ratio',ratio,...
    'flag_PANfromGT',flag_PANfromGT,'SNR_PANfromGT',snr_PANfromGT,...
    'request',{'MS_LR','PAN','GT','EXP'});

MS_LR=I_cell{1}; PAN=I_cell{2}; GT=I_cell{3}; EXP=I_cell{4}; clear I_cell;

compression_max=2^((ratio^2/(ratio^2+MS_LR.size(3)))*MS_LR.IntBits);

% Radiometric Resolution

PAN_radres=PAN;
MS_LR_radres=MS_LR;
PAN_radres.data=round(round(PAN.data/PAN.DynamicRange*compression_max)*PAN.DynamicRange/compression_max);
MS_LR_radres.data=round(round(MS_LR.data/MS_LR.DynamicRange*compression_max)*MS_LR.DynamicRange/compression_max);


MF_radres=Fusion_Algorithms_Pansharpening(...
    MS_LR_radres, PAN_radres,'GT',GT,'methods_list',methods_list);

[MR_radres,qindex_list,MR_idx_radres]=indexes_evaluation_RR(MF_radres,GT,'qindex_list', qindex_list);

% JPEG compression

PAN_jpeg=PAN;
MS_LR_jpeg=MS_LR;
PAN_jpeg.data=uint8(PAN.data/PAN.DynamicRange*255);
imwrite(PAN_jpeg.data,'temp.j2k','CompressionRatio',log2(255)/log2(compression_max));
PAN_jpeg.data=imread('temp.j2k');
PAN_jpeg.data=(round(double(PAN_jpeg.data)*PAN.DynamicRange/255));

MS_LR_jpeg.data=uint8(MS_LR.data/MS_LR.DynamicRange*255);

if MS_LR.size(3)==3
    imwrite(MS_LR_jpeg.data,'temp.j2k','CompressionRatio',log2(255)/log2(compression_max));
    MS_LR_jpeg.data=double(imread('temp.j2k'));
    MS_LR_jpeg.data=round(MS_LR_jpeg.data*MS_LR.DynamicRange/255);
else
    temp=[];
    for kk=1:size(MS_LR.data,3)
        imwrite(MS_LR_jpeg.data(:,:,kk),'temp.j2k','CompressionRatio',log2(255)/log2(compression_max));
        temp=cat(3,temp,double(imread('temp.j2k')));
    end
    MS_LR_jpeg.data=round(temp*MS_LR.DynamicRange/255);
end
    

MF_jpeg=Fusion_Algorithms_Pansharpening(...
    MS_LR_jpeg, PAN_jpeg,'GT',GT,'methods_list',methods_list);
[MR_jpeg,qindex_list,MR_idx_jpeg]=indexes_evaluation_RR(MF_jpeg,GT,'qindex_list', qindex_list);

MF=Fusion_Algorithms_Pansharpening(...
    MS_LR, PAN,'GT',GT,'methods_list',methods_list);
[MR,qindex_list,MR_idx_none]=indexes_evaluation_RR(MF,GT,'qindex_list', qindex_list);

mkdir(fullfile(output_folder,im_tag));

if flag_PANfromGT==1, string_sim=sprintf('sim%d',snr_PANfromGT); else, string_sim='real'; end
filename=sprintf([im_tag,'_r%d_',string_sim,'_compression'],ratio);
savefile=fullfile(output_folder,im_tag,filename);
save([savefile,'.mat'],'MF_radres','MR_radres','MF_jpeg','MR_jpeg','MF','MR','qindex_list')

qindex_metric=1;
idx_img=find(MR_idx_none(qindex_metric,:)==1);
idx_img_radres=find(MR_idx_radres(qindex_metric,:)==1);
idx_img_jpeg=find(MR_idx_jpeg(qindex_metric,:)==1);

a=viewimage(MF.data(:,:,MF.Bands_to_display,idx_img)); imshow(a,[],'Border','tight'); imwrite(a,[savefile,'_FUSION_BEST.png']);
a=viewimage(MF_radres.data(:,:,MF_radres.Bands_to_display,idx_img_radres)); imshow(a,[],'Border','tight'); imwrite(a,[savefile,'_COMPRADRES_BEST.png']);
a=viewimage(MF_jpeg.data(:,:,MF_jpeg.Bands_to_display,idx_img_jpeg)); imshow(a,[],'Border','tight'); imwrite(a,[savefile,'_COMPJPEG_BEST.png']);


MR_label=cell(1,3*numel(methods_list));
for ii=1:numel(methods_list)
    MR_label{ii}=methods_list{ii};
    MR_label{ii+numel(methods_list)}=['BIN+',methods_list{ii}];
    MR_label{ii+2*numel(methods_list)}=['JPEG+',methods_list{ii}];
end
MR_out=cat(2,MR,cat(2,MR_radres,MR_jpeg));
MR_idx=cat(2,MR_idx_none,cat(2,MR_idx_radres,MR_idx_jpeg));

matrix2latex(MR_out.','filename',[savefile,'.tex'],...
    'row',MR_label,'col',qindex_list,'align','c','significant',4,...
    'bold',MR_idx.'==1,'underline',MR_idx.'==2);

% I_vis=Visualization_image({GT,PAN,MS_LR,MF},'printMAT',1,'output_file','Test\output_test');

% string_results=matrix2latex_mod4(MR.','column',qindex_list,...
%    'row',methods_list,'flag_highlight',1,'output','Test\output_test','digits',4);

%I_out=Fusion_Algorithms_Pansharpening(I_MS_LR,I_PAN,'methods',fusion_list);


% % Radiometric Resolution
% 
% PAN_radres=PAN;
% MS_LR_radres=MS_LR;
% PAN_radres.data=round(round(PAN.data/PAN.DynamicRange*255)/255*PAN.DynamicRange);
% MS_LR_radres.data=round(round(MS_LR.data/MS_LR.DynamicRange*255)/255*PAN.DynamicRange);
% 
% 
% MF_radres=Fusion_Algorithms_Pansharpening(...
%     MS_LR_radres, PAN_radres,'GT',GT,'methods_list',methods_list);
% 
% [MR_radres,qindex_list]=indexes_evaluation_RR(MF_radres,GT,'qindex_list', qindex_list);
