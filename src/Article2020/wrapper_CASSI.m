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

function [I_out,I_acq,mask,MR]=wrapper_CASSI(varargin)

%% Support folders' paths
rng('default');  % For reproductible results
current_folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(current_folder,'..','..','data','output');
support_folder=fullfile(current_folder,'..','..','support');
addpath(fullfile(support_folder,'Load_info'),...
        fullfile(support_folder,'Interpolation'),...
        fullfile(support_folder,'Inversion_solver'),...
        fullfile(support_folder,'Mosaic'),...
        fullfile(support_folder,'Operator'),...
        fullfile(support_folder,'Validation'),...
        fullfile(support_folder,'Visualization'));

%% Default variables

im_tag='Washington_cut256_RGB';
testtype='default';
lambda_v=0.01;
preproc='regrnonneg';
ratio=[];
mask_tag='mindis';
alpha=[];
qindex_list={'SSIM','PSNR','ERGAS','SAM','sCC','UIQI','Q2n'};
idx_metric=1;
inversion='TV';
Nbiter=1000;
tol=0;
flag_vis=1;
flag_PANfromGT=0;
SNR_PANfromGT=[];
ra=[];
output_fol = 'test_compression';


%% Variables setter

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1}; 
    if any(strcmpi(pname,{'im_tag','image','img','im'}))
        im_tag=pval;
    elseif any(strcmpi(pname,{'mask','mask_tag'}))          % Tag for the used mask type
        mask_tag=pval;
    elseif any(strcmpi(pname,{'alpha'}))                    % Percentage of the mask assigned to the HR image
        alpha=pval;
    elseif any(strcmpi(pname,{'lambda','lambda_v'}))        % vector of the regularization parameters
        lambda_v=pval;
    elseif any(strcmpi(pname,{'testtype','test'}))          % Id of the type of test to perform
        testtype=pval;
    elseif any(strcmpi(pname,{'preproc','preprocessing'}))  % Preprocessing algorithm
        preproc=pval;
    elseif any(strcmpi(pname,{'ratio','scale'}))            % Scale ratio between PAN and MS
        ratio=pval;
    elseif any(strcmpi(pname,{'qindex_list','qindex','index','index_list'})) % List of quality indices for validation
        qindex_list=pval;
    elseif any(strcmpi(pname,{'idx_metric','metric'}))      % Index of the metric to use for best lambda
        idx_metric=pval;
    elseif any(strcmpi(pname,{'inv','inversion'}))          % Inversion method tag
        inversion=pval;
    elseif any(strcmpi(pname,{'Nbiter','iter'}))            % Number of iterations for the recursive inversion algorithm
        Nbiter=pval;
    elseif any(strcmpi(pname,{'tol','tolerance'}))          % Tolerance for the error in the iteration
        tol=pval;
    elseif any(strcmpi(pname,{'figure','visualization','vis'}))   % Flag for visualizing images
        flag_vis=pval;
    elseif any(strcmpi(pname,{'flag_PANfromGT','simulated','sim','PANfromGT'}))
        flag_PANfromGT=pval;
    elseif any(strcmpi(pname,{'SNR_PANfromGT','SNR'}))
        SNR_PANfromGT=pval; 
    elseif any(strcmpi(pname,{'radius','ra','magnify'}))
        ra=pval;
    elseif any(strcmpi(pname,{'output','folder','output_folder','output_fol'}))
        output_fol=pval;
    end
end

output_folder = fullfile(output_folder, output_fol);

%% Load Methods
if strcmpi(testtype,'msonly')
    ratio_mask=1;
    ratio=1;
    mask_tag_PAN='none';
elseif strcmpi(testtype,'nodegrad')
    error('Modality is currently not supported');
    % ratio_mask=ratio; 
    % ratio=1;
    % mask_tag_PAN=[];
elseif strcmpi(testtype,'nomask')
    ratio_mask=ratio;
    mask_tag='none';
    mask_tag_PAN='none';
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
    I_PAN.spectralweights=[];
    I_PAN.DynamicRange=I_EXP.DynamicRange;
end

%% Setting Operator for masking operation

mask=load_mask('sizes',I_EXP.size,'type',mask_tag,'mask_PAN',mask_tag_PAN,'band_HR',I_PAN.size(3),'ratio',ratio_mask,'start_pos',I_EXP.start_pos,'alpha',alpha);
im.HR=I_PAN.data; im.LR=I_EXP.data;
y=mask.mosaic(im,mask); % Acquisition (NOTE: No pre-processing at all!)

% opA =load_degmaskoperator('lpfilter',I_EXP.KerBlu,'spectralweights',I_EXP.spectralweights,'mask',mask,'flag_fft',0);    
% y = opA.mosaic(cat(3,I_PAN.data,I_EXP.data)); % Acquisition (NOTE: No pre-processing at all!)

%% Estimation of the HR image degradation model
I_demo = mask.demosaic(y,mask);

I_EXP_WB=interp2D_mosaic(I_demo.image_EXP,'mask',I_demo.mask_EXP,'method','WB','nn',100);
%I_EXP_demo  = demosaic_classic(I_demo.image_EXP,I_demo.mask_EXP,'WB'); I_EXP_WB=I_EXP_demo.data;

if strcmpi(testtype,'msonly')
    I_PAN_WB=[];
    I_PAN_LR=[];
    I_PAN_temp=permute(sum(I_EXP_WB.*shiftdim(I_EXP.spectralweights,-2),3),[1,2,4,3]);
    std_PAN=std(I_PAN_temp,0,1:3);
    mean_PAN=mean(I_PAN_temp,1:2);
else
    I_PAN_WB=interp2D_mosaic(I_demo.image_HR,'mask',I_demo.mask_HR,'method','WB','nn',100);
    % I_PAN_demo=demosaic_classic(I_demo.mosaic_HR,I_demo.mask_HR,'WB'); I_PAN_WB=I_PAN_demo.data;
    I_PAN_LR=imresize(imresize(I_PAN_WB,1/ratio),ratio);
    std_PAN=std(I_PAN_WB,0,1:3);
    mean_PAN=mean(I_PAN_WB,1:2);
end

if strncmpi(preproc,'regr',4)
    I_EXP_norm=(I_EXP_WB-mean(I_EXP_WB,1:2))./std(I_EXP_WB,0,1:2).*std_PAN+mean_PAN;
    I_PAN_norm=I_PAN_WB;
    weights=zeros(size(I_EXP_WB,3),size(I_PAN_WB,3));
    for kk=1:size(I_PAN_WB,3)  
        if strcmpi(preproc,'regr'), weights(:,kk)=estimation_alpha(I_EXP_WB,I_PAN_LR(:,:,kk),'global'); end
        if strcmpi(preproc,'regrnonneg'), weights(:,kk)=lsqnonneg(reshape(I_EXP_WB,[],size(I_EXP_WB,3)),reshape(I_PAN_LR(:,:,kk),[],1)); end
        if strcmpi(preproc,'regrsum1')
            weights(:,kk) = lsqlin(reshape(I_EXP_WB,[],size(I_EXP_WB,3)),...
                reshape(I_PAN_LR(:,:,kk),[],1),[],[],ones(size(I_EXP_WB,3)),1,0,[]);
        end
        if strcmpi(preproc,'regravg'), weights(:,kk)=1/size(I_EXP_WB,3); end
    end
    max_value=I_PAN.DynamicRange;
elseif strcmpi(preproc,'hism')
    I_EXP_norm=(I_EXP_WB-mean(I_EXP_WB,1:2))./std(I_EXP_WB,0,1:2).*std_PAN+mean_PAN;
    I_PAN_norm=I_PAN_WB;
    weights=I_PAN.spectralweights;
    max_value=I_PAN.DynamicRange;
elseif strcmpi(preproc,'none')
    I_EXP_norm=I_EXP_WB;
    I_PAN_norm=I_PAN_WB;
    weights=I_PAN.spectralweights;
    max_value=I_EXP.DynamicRange;
end

fprintf('Regression Weights:');
for ii=1:length(weights), fprintf(' %.4f',weights(ii)); end; fprintf('\n');

% I_EXP_norm=(I_EXP_WB-subtract_MS)./divide_MS;
% I_PAN_norm=(I_PAN_WB-subtract_PAN)./divide_PAN;
opA =load_degmaskoperator('lpfilter',I_EXP.KerBlu,'spectralweights',weights,...
    'mask',mask,'flag_fft',0,'radius',ra);

y=opA.mosaic(cat(3,I_PAN_norm,I_EXP_norm));

%% Other operators
if iscell(inversion)
    opL_tag=inversion{1};
    opg_tag=inversion{2};
    opW_tag=inversion{3};
elseif strcmpi(inversion,'TV')
    opL_tag='TVc';
    opg_tag='Norm_l221';
    opW_tag='none';
elseif strcmpi(inversion,'Wavelet')
    opL_tag='none';
    opg_tag='Norm_l111';
    opW_tag='CAS8_sym8';
end

opg=load_operator(opg_tag);
opL=load_operator(opL_tag);
opW=load_operator(opW_tag,size(opA.adjoint(y_norm)));


%% Test for Regularization Parameter (lambda) 
% Note: Algorithm converges very slowly with very low or absent 
% levels of additive Gaussian noise

tau=1;
rho=1.9;

[x_lambda,cost_temp]=Loris(y_norm,'opA',opA,'opL',opL,'opW',opW,...
    'opg',opg,'rho',rho,'tau',tau,'lambda',lambda_v*max_value,...
    'init',opA.adjoint(y_norm),'Nbiter',Nbiter,'tol',tol);

% Recover moments
if ~any(strcmpi(preproc,{'none'}))
    x_lambda=(x_lambda-mean_PAN)./std_PAN.*std(I_EXP_WB,0,1:2)+mean(I_EXP_WB,1:2);
    % x_lambda=(x_lambda.*divide_MS)+subtract_MS;
end

fprintf('Quality index calculation:\n')
[MR_out,qi_label,MR_idx]= indexes_evaluation_RR(x_lambda,I_GT,'qindex_list',qindex_list);

if flag_PANfromGT==1, sim_string=sprintf('sim%d',round(SNR_PANfromGT)); else, sim_string='real'; end
if isempty(ra), ra=1; end


output_folder=fullfile(output_folder,im_tag);
mkdir(output_folder);
filename=sprintf('%s_r%d_%s_%s_%s_%s_m%d_%s_%s',im_tag,ratio,sim_string,mask_tag,testtype,preproc,round((ra-1)*10),opL_tag,opg_tag);

MR_label=cell(1,size(MR_out,2)); for ii=1:size(MR_out,2), MR_label{ii}=sprintf('lambda= %.2E',lambda_v(ii)); end
matrix2latex(MR_out.','filename',fullfile(output_folder,[filename,'.tex']),...
    'row',MR_label,'col',qi_label,'align','c','significant',4,...
    'bold',MR_idx.'==1,'underline',MR_idx.'==2);

[~,idx_lambda]=min(MR_idx(idx_metric,:));
fprintf('Mean of GT: %.4f, Mean of Result: %.4f\n',mean(I_EXP.data(:)),mean(reshape(x_lambda(:,:,:,idx_lambda),numel(I_GT.data),[])));

%% Visualization

if flag_vis==1
    fig_idx=1; fig_idx2=2;
    figure(fig_idx2); imshow(1,[]); hold off;
    if ~strcmpi(testtype,'nomask')
        a=visualize_mask(mask,0); figure(fig_idx2); imshow(a,[],'Border','tight'); title('MASK'); imwrite(a,fullfile(output_folder,[filename,'_MASK.png']));
        figure(fig_idx); subplot(244); a=viewimage_outputonly(y);  imshow(a,[],'Border','tight'); title('Acquisition'); imwrite(a,fullfile(output_folder,[filename,'_COMP.png']));
    end
    if ~strcmpi(testtype,'msonly')
        figure(fig_idx); subplot(241); a=viewimage_outputonly(I_PAN.data(:,:,I_PAN.Bands_to_display)); imshow(a,[],'Border','tight'); title('PAN'); imwrite(a,fullfile(output_folder,[filename,'_PAN.png']));
    end
    figure(fig_idx); subplot(242); a=viewimage_outputonly(I_EXP.data(:,:,I_EXP.Bands_to_display)); imshow(a,[],'Border','tight'); title('EXP'); imwrite(a,fullfile(output_folder,[filename,'_EXP.png']));
    figure(fig_idx); subplot(243); a=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display)); imshow(a,[],'Border','tight'); title('GT'); imwrite(a,fullfile(output_folder,[filename,'_GT.png']));
    init=opA.adjoint(y);
    figure(fig_idx); subplot(245); a=viewimage_outputonly(init(:,:,I_MS_LR.Bands_to_display)); imshow(a,[],'Border','tight'); title('Initialization'); imwrite(a,fullfile(output_folder,[filename,'_INIT.png']));
    figure(fig_idx); subplot(246); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,1)); imshow(a,[],'Border','tight'); title(sprintf('Minimum (\\lambda=%.2E)',lambda_v(1))); imwrite(a,fullfile(output_folder,[filename,'_INVMIN.png']));
    figure(fig_idx); subplot(247); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,idx_lambda)); imshow(a,[],'Border','tight'); title(sprintf('Best (\\lambda=%.2E)',lambda_v(idx_lambda))); imwrite(a,fullfile(output_folder,[filename,'_INVBEST.png']));
    figure(fig_idx); subplot(248); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,end)); imshow(a,[],'Border','tight'); title(sprintf('Maximum (\\lambda=%.2E)',lambda_v(end))); imwrite(a,fullfile(output_folder,[filename,'_INVMAX.png']));
    
    figure; hold off; title(sprintf('%s Evaluation',qi_label{idx_metric}));
    ylabel('Value'); xlabel('Regularization parameter (\lambda)');
    plot(lambda_v,MR_out(idx_metric,:));
    drawnow;
end


%% Generation of outputs

I_out=I_GT;
I_out.data=x_lambda;
I_out.test=testtype;
I_out.inversion={opL_tag,opg_tag,opW_tag};
I_out.preprocessing=preproc;
I_out.label='Fused';
I_out.spectralweights=I_EXP.spectralweights;
I_out.KerBlu=I_EXP.KerBlu;
I_out.lambda=lambda_v;
I_out.iterations=Nbiter;
I_out.mask=mask_tag;
I_out.maskalpha=alpha;
I_out.cost=cost_temp;
I_out.magnifyradius=ra;

I_acq=I_GT;
I_acq.data=y;
I_acq.type='PAN';
I_acq.label='Acquisition';
I_acq.test=testtype;
I_acq.spectralweights=I_EXP.spectralweights;
I_acq.KerBlu=I_EXP.KerBlu;
I_acq.preproc=preproc;

MR.data=MR_out;
MR.qindex=qi_label;
MR.label=MR_label;
MR.bestindex=MR_idx;

save(fullfile(output_folder,[filename,'.mat']),'I_out','I_acq','MR','mask','I_load');