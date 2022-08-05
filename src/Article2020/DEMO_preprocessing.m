%% Adding support path

clearvars; close all;
rng('default');  % For reproductible results
current_folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(current_folder,'..','..','data','output','test_compression');
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'load_info'));
addpath(fullfile(project_folder,'operator'));
addpath(fullfile(project_folder,'inversion_solver'));
addpath(fullfile(project_folder,'mosaic'));
addpath(fullfile(project_folder,'visualization'));

%% Parameters
% preproc='none'; 
% preproc='hism'; 
% preproc='regr'; 
% preproc='regrnobias';
% preproc='regrnonneg';
im_tag='Washington_cut256_RGB';
lambda_v=logspace(-3,-1,5);
testtype='TV';
ratio=2;
mask_tag='mindis'; mask_tag_PAN='default';
alpha=0.5;
qindex_list={'psnr','ERGAS','SAM','sCC','UIQI','q2n','SSIM'}; idx_metric=7;

%% Load input image
I_load=Load_Dataset_Pansharpening(im_tag,'request',{'MS_LR','PAN','EXP','GT'},'ratio',ratio,'flag_PANfromGT',0);
I_MS_LR=I_load{1}; I_PAN=I_load{2}; I_EXP=I_load{3}; I_GT=I_load{4};

%% Generation of simulated CASSI acquisition
if strcmpi(preproc,'norm')
    im=cat(3,I_PAN.data/I_PAN.DynamicRange,I_EXP.data/I_EXP.DynamicRange);
    max_value=1; % max_value=I_load{3}.DynamicRange; % If not normalized
elseif strcmpi(preproc,'hism')
    I_PAN_LR=imresize(I_PAN.data,1/ratio);
    I_EXP.data=(I_EXP.data-mean(I_EXP.data,1:2))./std(I_EXP.data,0,1:2).*std(I_PAN.data,0,1:2)+mean(I_PAN.data,1:2);
    im=cat(3,I_PAN.data/I_PAN.DynamicRange,I_EXP.data/I_EXP.DynamicRange);
    max_value=1; % max_value=I_load{3}.DynamicRange; % If not normalized
elseif strcmpi(preproc,'regr')
    addpath(fullfile(project_folder,'filter'));
    I_PAN_LR=imresize(imresize(I_PAN.data,1/ratio),ratio);
    weights=estimation_alpha2(I_EXP.data,I_PAN_LR,'global');
    fprintf('Regression Weights:');for ii=1:length(weights), fprintf(' %.4f',weights(ii)); end; fprintf('\n');
    I_EXP.spectralweights=weights(1:end-1);
    I_PAN.data=I_PAN.data-weights(end);
    im=cat(3,I_PAN.data/I_PAN.DynamicRange,I_EXP.data/I_EXP.DynamicRange);
    max_value=1; % max_value=I_load{3}.DynamicRange; % If not normalized
elseif strcmpi(preproc,'regrnobias')
    addpath(fullfile(project_folder,'filter'));
    I_PAN_LR=imresize(imresize(I_PAN.data,1/ratio),ratio);
    weights=estimation_alpha(I_EXP.data,I_PAN_LR,'global');
    fprintf('Regression Weights:');for ii=1:length(weights), fprintf(' %.4f',weights(ii)); end; fprintf('\n');
    I_EXP.spectralweights=weights;
    im=cat(3,I_PAN.data/I_PAN.DynamicRange,I_EXP.data/I_EXP.DynamicRange);
    max_value=1; % max_value=I_load{3}.DynamicRange; % If not normalized
elseif strcmpi(preproc,'regrnonneg')
    addpath(fullfile(project_folder,'filter'));
    I_PAN_LR=imresize(imresize(I_PAN.data,1/ratio),ratio);
    weights = lsqnonneg(reshape(I_EXP.data,[],I_EXP.size(3)),I_PAN_LR(:));
    fprintf('Regression Weights:');for ii=1:length(weights), fprintf(' %.4f',weights(ii)); end; fprintf('\n');
    I_EXP.spectralweights=weights;
    im=cat(3,I_PAN.data/I_PAN.DynamicRange,I_EXP.data/I_EXP.DynamicRange);
    max_value=1; % max_value=I_load{3}.DynamicRange; % If not normalized
elseif strcmpi(preproc,'none')
    im=cat(3,I_PAN.data,I_EXP.data);
    max_value=I_EXP.DynamicRange; % max_value=I_load{3}.DynamicRange; % If not normalized
end


%% Setting Operator for masking operation

mask=load_mask('sizes',I_GT.size,'mask_MS',mask_tag,'mask_PAN',mask_tag_PAN,'band_HR',I_PAN.size(3),'ratio',ratio,'alpha',alpha);

opA =load_degmaskoperator('lpfilter',I_EXP.KerBlu,'spectralweights',I_EXP.spectralweights,'mask',mask,'flag_fft',0);
opg=load_operator('Norm_l221');
opL=load_operator('TVc');
opW=load_operator('none');

y = opA.masking(im); % Acquisition


%% Test for Regularization Parameter (lambda) 
% Note: Algorithm converges very slowly with very low or absent 
% levels of additive Gaussian noise

% lambda_v=logspace(-3,-1,3);
% lambda_v=0.001;
tau=1;
rho=1.9;
Nbiter=1000;
tol=0;

[x_lambda,cost_temp]=Loris(y,'opA',opA,'opL',opL,'opW',opW,'opg',opg,...
    'rho',rho,'tau',tau,'lambda',lambda_v*max_value,'init',opA.adjoint(y),'Nbiter',Nbiter,...
    'tol',tol);

if any(strcmpi(preproc,{'norm','regr','regrnobias','regrnonneg'}))    
    x_lambda=x_lambda*I_EXP.DynamicRange;
elseif strcmpi(preproc,'hism')
    x_lambda=x_lambda*I_EXP.DynamicRange;
    x_lambda=(x_lambda-mean(I_PAN.data,1:2))./std(I_PAN.data,0,1:2).*std(I_EXP.data,0,1:2)+mean(I_EXP.data,1:2);
end

fprintf('Quality index calculation:\n')
[MR,qi_label]= indexes_evaluation_RR(x_lambda,I_GT,'qindex_list',qindex_list);


figure; hold on; title('SSIM Evaluation');
ylabel('Value'); xlabel('Regularization parameter (\lambda)');
plot(lambda_v,MR(idx_metric,:));

[~,idx_lambda]=max(MR(idx_metric,:));
fprintf('Mean GT (before): %.4f, (after): %.4f\n',mean(I_EXP.data(:)),mean(reshape(x_lambda(:,:,:,idx_lambda),numel(I_GT.data),[])));

filename=sprintf('%s_%s_%s_%s',im_tag,mask_tag,preproc,testtype);
save(fullfile(output_folder,[filename,'.mat']));

fig_idx=10; fig_idx2=11;
figure(fig_idx2+1); a=visualize_mask(mask,0); imshow(a,[],'Border','tight');
figure(fig_idx); subplot(241); a=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display)); imshow(a,[],'Border','tight'); title('GT'); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_GT']),'png');
figure(fig_idx); subplot(242); a=viewimage_outputonly(I_EXP.data(:,:,I_EXP.Bands_to_display)); imshow(a,[],'Border','tight'); title('EXP'); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_EXP']),'png');
figure(fig_idx); subplot(243); a=viewimage_outputonly(I_PAN.data(:,:,I_PAN.Bands_to_display)); imshow(a,[],'Border','tight'); title('PAN'); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_PAN']),'png');
figure(fig_idx); subplot(244); a=viewimage_outputonly(y);  imshow(a,[],'Border','tight'); title('Acquisition'); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_COMP']),'png');
init=opA.adjoint(y);
figure(fig_idx); subplot(245); a=viewimage_outputonly(init(:,:,I_MS_LR.Bands_to_display)); imshow(a,[],'Border','tight'); title('Initialization'); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_INIT']),'png');
figure(fig_idx); subplot(246); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,1)); imshow(a,[],'Border','tight'); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(1))); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_INVMIN']),'png');
figure(fig_idx); subplot(247); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,idx_lambda)); imshow(a,[],'Border','tight'); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(idx_lambda))); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_INVBEST']),'png');
figure(fig_idx); subplot(248); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,end)); imshow(a,[],'Border','tight'); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(end))); figure(fig_idx2); imshow(a,[],'Border','tight'); saveas(gcf,fullfile(output_folder,[filename,'_INVMAX']),'png');
