current_folder=fileparts(mfilename('fullpath'));
support_folder=fullfile(current_folder,'..','support');
addpath(fullfile(support_folder,'load_info'),...
        fullfile(support_folder,'fusion'),...
        fullfile(support_folder,'visualization'));

% im_tag='Hobart1_RGB';
% im_tag='Hobart2_RGB';
% im_tag='Rome_RGB';
% im_tag='Sofia_cut2_ALI_ALI_RGB';
% im_tag='Hobart2';
% im_tag='China';
% im_tag='RdJ_WV3_WV3_4';
im_tag='RdJ_WV3_WV3';
% im_tag='Rio_WV2_WV2';
% im_tag='Rome_WV2_WV2';

ratio=2;
flag_PANfromGT=0;
methods_list={'EXP','GSA','MTF-GLP-CBD','BDSD','BayesNaive'};
qindex_list={'SSIM','PSNR','Q2^n','SAM','ERGAS','SCC','UIQI'};


I_cell=Load_Dataset_Pansharpening_RR(im_tag,'ratio',ratio,...
    'flag_PANfromGT',flag_PANfromGT,'request',{'MS_LR','PAN','GT','EXP'});

% MS_LR=I_cell{1}; PAN=I_cell{2}; GT=I_cell{3}; EXP=I_cell{4}; clear I_cell;

% MF=Fusion_Algorithms_Pansharpening(...
%     MS_LR,PAN,'EXP',EXP,'GT',GT,'methods_list',methods_list);

% [MR,qindex_list]=indexes_evaluation_RR(MF,GT,'qindex_list', qindex_list);

Visualization_image(I_cell,'flag_cutbounds',1,'dim_cut',5,'flag_thvalues',1,'printMAT',1,'output_file','Test\output_test');



