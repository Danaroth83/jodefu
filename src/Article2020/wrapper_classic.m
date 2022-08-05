%% COMPRESSED ACQUISITIONS
%
% Author:
% Daniele Picone
%
% Description:
% This wrapper allows to simulate a compressed acquisition with a digital 
% micromirror device (DMD) which mosaicks acquisitions from both a high
% resolution (HR) monoband image and a low resolution (LR) multiband image
% and allows to invert the acquisition to generate a fused image via a
% variational algorithm
%
% Usage:
% [I_out,I_acq,mask,MR]=wrapper_compressedacquisition('Field',Value)
%
% Output:
% I_out: Struct of the result of the fusion via variational algorithm
% I_acq: Struct of the simulated acquisition
% mask:  Struct of the mask implemented by the DMD
% MR:    Struct of the results of the Reduced Resolution tests
%
% Fields:
% 'im': tag for the set of HR and LR image to open
%       (A list of tags is available in ../../support/Load_info/Tags_list.txt)
% 'sim': If 1, simulates the HR image from the ground truth (default: 0)
% 'ratio': Scale ratio between the HR and LR image (default: Ratio of the original images)
% 'qindex': cell of strings describing the quality indices
%           (default: {'SSIM','PSNR','ERGAS','SAM','sCC','UIQI','Q2n'})
% 'metric': Index of the cell above to decide the best lambda (default: 1)
% 'vis': If 1, shows images on screen and prints figures (default: 1)

function [I_out,I_acq,mask,MR]=wrapper_classic(varargin)

%% Support folders' paths
rng('default');  % For reproductible results
current_folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(current_folder,'..','..','data','output');
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'load_info'),...
        fullfile(project_folder,'demosaic'),...
        fullfile(project_folder,'fusion'),...
        fullfile(project_folder,'interpolation'),...
        fullfile(project_folder,'mosaic'),...
        fullfile(project_folder,'validation'),...
        fullfile(project_folder,'visualization'));

%% Default variables

im_tag='Washington_cut256_RGB';
ratio=[];
mask_tag='Bayer';
alpha=[];
qindex_list={'SSIM','PSNR','ERGAS','SAM','sCC','UIQI','Q2n'};
idx_metric=1;
flag_vis=1;
flag_PANfromGT=0;
SNR_PANfromGT=[];
interpolation='WB';
demosaic_list='ID';
fusion_list='GSA';
testtype='default';
output_fol='test_compression';

%% Variables setter

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1}; 
    if any(strcmpi(pname,{'im_tag','image','img','im'}))
        im_tag=pval;
    elseif any(strcmpi(pname,{'mask','mask_tag'}))          % Tag for the used mask type
        mask_tag=pval;
    elseif any(strcmpi(pname,{'testtype','test'}))          % Id of the type of test to perform
        testtype=pval;
    elseif any(strcmpi(pname,{'alpha'}))                    % Percentage of the mask assigned to the HR image
        alpha=pval;
    elseif any(strcmpi(pname,{'interpolation','method_interpolation'}))
        interpolation=pval;
    elseif any(strcmpi(pname,{'demosaic','method_demosaic'}))
        demosaic_list=pval;
    elseif any(strcmpi(pname,{'fusion','method_fusion'}))
        fusion_list=pval;
    elseif any(strcmpi(pname,{'ratio','scale'}))            % Scale ratio between PAN and MS
         ratio=pval;
    elseif any(strcmpi(pname,{'qindex_list','qindex','index','index_list'})) % List of quality indices for validation
        qindex_list=pval;
    elseif any(strcmpi(pname,{'idx_metric','metric'}))      % Index of the metric to use for best lambda
        idx_metric=pval;
    elseif any(strcmpi(pname,{'figure','visualization','vis'}))   % Flag for visualizing images
        flag_vis=pval;
    elseif any(strcmpi(pname,{'flag_PANfromGT','simulated','sim','PANfromGT'}))
        flag_PANfromGT=pval;
    elseif any(strcmpi(pname,{'SNR_PANfromGT','SNR'}))
        SNR_PANfromGT=pval;
    elseif any(strcmpi(pname,{'output','folder','output_folder','output_fol'}))
        output_fol=pval;
    end
end

output_folder = fullfile(output_folder, output_fol);

if iscell(interpolation),  interpolation=interpolation{1}; end
if ~iscell(demosaic_list), demosaic_list={demosaic_list}; end
if ~iscell(fusion_list),   fusion_list={fusion_list}; end

%% Load Methods

if strcmpi(testtype,'msonly')
    ratio_mask=1;
    ratio=1;
    mask_tag_PAN='none';
    fusion_list={'none'};
elseif strcmpi(testtype,'nodegrad')
    error('Modality is currently not supported');
    % ratio_mask=ratio; 
    % ratio=1;
    % mask_tag_PAN=[];
    % fusion_list={'none'};
elseif strcmpi(testtype,'nomask')
    ratio_mask=ratio;
    mask_tag='none';
    mask_tag_PAN='none';
    demosaic_list={'none'};
else
    ratio_mask=ratio;
    mask_tag_PAN=[];
end

%% Load input image
I_load=Load_Dataset_Pansharpening(im_tag,'request',{'MS_LR','PAN','EXP','GT'},...
'ratio',ratio,'flag_PANfromGT',flag_PANfromGT,'SNR_PANfromGT',SNR_PANfromGT);
I_MS_LR=I_load{1}; I_PAN=I_load{2}; I_EXP=I_load{3}; I_GT=I_load{4};
ratio=I_PAN.ratio;

%% Fix to impose zero HR bands
if strcmpi(testtype,'msonly') % Kinda hacky, but works
    I_PAN.data=[];
    I_PAN.size=[0,0,0];
end

%% Setting Operator for masking operation

mask=load_mask('sizes',I_EXP.size,'type',mask_tag,'mask_PAN',mask_tag_PAN,...
    'band_HR',I_PAN.size(3),'ratio',ratio_mask,'start_pos',I_EXP.start_pos,'alpha',alpha);
im.HR=I_PAN.data; im.LR=I_EXP.data;
y=mask.mosaic(im,mask); % Acquisition (NOTE: No pre-processing at all!)


%% Demosaic + Fusion Procedure
I_demo=mask.demosaic(y,mask);

% PAN inpainting
% I_PAN_demo  = interp2D_mosaic(I_demo.image_HR,'mask',I_demo.mask_HR,'method',interpolation,'nn',100);
I_PAN_demo  = demosaic_classic(I_demo.mosaic_HR,I_demo.mask_HR,interpolation); % Some extra fields are not necessary

% MS demosaic
I_MS_demo  = demosaic_classic(I_demo.mosaic_LR,I_demo.mask_LR,demosaic_list,I_MS_LR.DynamicRange,I_MS_LR.wavelength);
[MR_MS,q_index_MS,MR_MS_idx]= indexes_evaluation_RR(I_MS_demo.data,I_MS_LR.data,'qindex_list',qindex_list,'ratio',ratio);
[~,idx_best_MS]=min(MR_MS_idx(idx_metric,:));
fprintf('Best %s: %.4f. Method: %s\n',q_index_MS{idx_metric}, MR_MS(idx_metric,idx_best_MS),demosaic_list{idx_best_MS});
    
% Fusion

x=zeros(size(I_MS_demo.data,1)*ratio,size(I_MS_demo.data,2)*ratio,size(I_MS_demo.data,3),numel(demosaic_list),numel(fusion_list));
for kk=1:numel(demosaic_list)
    I_MS_in=I_MS_LR;
    I_MS_in.data=I_MS_demo.data(:,:,:,kk);
    I_PAN_in=I_PAN;
    I_PAN_in.data=I_PAN_demo.data;
    I_out_temp=Fusion_Algorithms_Pansharpening(I_MS_in,I_PAN_in,'methods',fusion_list);
    x(:,:,:,kk,:)=I_out_temp.data;
end
x=reshape(x,size(x,1),size(x,2),size(x,3),[]);

%% Validation
fprintf('Quality index calculation:\n')
[MR_out,qi_label,MR_idx]= indexes_evaluation_RR(x,I_GT,'qindex_list',qindex_list);

MR_label=strings(1,size(MR_out,2)*size(MR_out,3)); 
for jj=1:numel(fusion_list)
    for ii=1:numel(demosaic_list)
        idx=ii+(jj-1)*numel(demosaic_list);
        MR_label{idx}=sprintf('%s/%s',demosaic_list{ii},fusion_list{jj});
    end
end

MR_out=reshape(MR_out,size(MR_out,1),[]);
MR_idx=reshape(MR_idx,size(MR_idx,1),[]);

if flag_PANfromGT==1, sim=sprintf('sim%d',round(SNR_PANfromGT)); else, sim='real'; end

output_folder=fullfile(output_folder,im_tag);
mkdir(output_folder);
filename=sprintf('%s_r%d_%s_%s_%s_fusdem_%s',im_tag,ratio,sim,mask_tag,testtype,interpolation);

matrix2latex(MR_out.','filename',fullfile(output_folder,[filename,'.tex']),...
'row',MR_label,'col',qi_label,'align','c','significant',4,...
'bold',MR_idx.'==1,'underline',MR_idx.'==2); 

idx_best=MR_idx(idx_metric,:)==1;
idx_seco=MR_idx(idx_metric,:)==2;
fprintf('Best %s: %.4f\n',qi_label{idx_metric},MR_out(idx_metric,idx_best));
fprintf('Mean of GT: %.4f, Mean of Result: %.4f\n',mean(I_EXP.data(:)),mean(reshape(x(:,:,:,idx_best),numel(I_GT.data),[])));

%% Visualization

if flag_vis==1
    fig_idx=1; fig_idx2=2;
    if ~strcmpi(testtype,'nomask')
        figure(fig_idx2); a=visualize_mask(mask,0); imshow(a,[],'Border','tight'); imwrite(a,fullfile(output_folder,[filename,'_MASK.png']));
    end
    figure(fig_idx); subplot(243); a=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display)); imshow(a,[],'Border','tight'); title('GT'); imwrite(a,fullfile(output_folder,[filename,'_GT.png']));
    figure(fig_idx); subplot(241); a=viewimage_outputonly(I_EXP.data(:,:,I_EXP.Bands_to_display)); imshow(a,[],'Border','tight'); title('EXP'); imwrite(a,fullfile(output_folder,[filename,'_EXP.png']));
    if ~strcmpi(testtype,'msonly')
        figure(fig_idx); subplot(242); a=viewimage_outputonly(I_PAN.data(:,:,I_PAN.Bands_to_display)); imshow(a,[],'Border','tight'); title('PAN'); imwrite(a,fullfile(output_folder,[filename,'_PAN.png']));
        figure(fig_idx); subplot(246); a=viewimage_outputonly(I_PAN_demo.data(:,:,I_PAN.Bands_to_display));imshow(a,[],'Border','tight'); title('PAN Demosaic'); imwrite(a,fullfile(output_folder,[filename,'_PANDEMO.png']));
    end
    if ~strcmpi(testtype,'nomask')
        figure(fig_idx); subplot(244); a=viewimage_outputonly(y); imshow(a,[],'Border','tight'); title('Acquisition'); imwrite(a,fullfile(output_folder,[filename,'_COMP.png']));
    end
    figure(fig_idx); subplot(245); a=viewimage_outputonly(imresize(I_MS_demo.data(:,:,I_EXP.Bands_to_display,idx_best_MS),[size(I_MS_demo.data,1),size(I_MS_demo.data,2)].*ratio)); imshow(a,[],'Border','tight'); title('MS Demosaic'); imwrite(a,fullfile(output_folder,[filename,'_MSDEMO.png']));
    figure(fig_idx); subplot(248); a=viewimage_outputonly(x(:,:,I_MS_LR.Bands_to_display,idx_seco)); imshow(a,[],'Border','tight'); title([MR_label{idx_seco},' (2nd)']); imwrite(a,fullfile(output_folder,[filename,'_INV2ND.png']));
    figure(fig_idx); subplot(247); a=viewimage_outputonly(x(:,:,I_MS_LR.Bands_to_display,idx_best)); imshow(a,[],'Border','tight'); title([MR_label{idx_best},' (Best)']); imwrite(a,fullfile(output_folder,[filename,'_INVBEST.png']));
    drawnow;
end


%% Generation of outputs

I_out=I_GT;
I_out.data=x;
I_out.interpol_HR=interpolation;
I_out.methods_demosaic=demosaic_list;
I_out.methods_fusion=fusion_list;
I_out.label='Fused';
I_out.test=testtype;
I_out.spectralweights=I_EXP.spectralweights;
I_out.KerBlu=I_EXP.KerBlu;
I_out.mask=mask_tag;
I_out.maskalpha=alpha;

I_acq=I_GT;
I_acq.data=y;
I_acq.type='PAN';
I_acq.label='Acquisition';
I_acq.test=testtype;
I_acq.spectralweights=I_EXP.spectralweights;
I_acq.KerBlu=I_EXP.KerBlu;

MR.data=MR_out;
MR.qindex=qi_label;
MR.label=MR_label;
MR.bestindex=MR_idx;

save(fullfile(output_folder,[filename,'.mat']),'I_out','I_acq','MR','mask','I_load');