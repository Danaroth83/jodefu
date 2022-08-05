%% Adding support path

clearvars; close all;
rng('default');  % For reproductible results
support_folder=fullfile('..','..','support');
output_folder=fullfile('..','..','data','output','test_compression');
addpath(fullfile(support_folder,'Load_info'));
addpath(fullfile(support_folder,'Operator'));
addpath(fullfile(support_folder,'Quality_indices'));
addpath(fullfile(support_folder,'Inversion_solver'));
addpath(fullfile(support_folder,'Mosaic'));
addpath(fullfile(support_folder,'Visualization'));


%% Loading Image and Generating Acquisition
im_tag='Washington_cut256_RGB';
mask_tag='CASSI';
ratio=2;

%% Loading Image and Generating Acquisition
I_load=Load_Dataset_Pansharpening(im_tag,'request',{'MS_LR','PAN','EXP','GT'},'ratio',ratio);
I_MS_LR=I_load{1}; I_PAN=I_load{2}; I_EXP=I_load{3}; I_GT=I_load{4};

%% Setting Operator for masking operation
mask=load_mask('sizes',I_EXP.size,'mask','CASSI','band_HR',I_PAN.size(3),'ratio',ratio,'alpha',0.5);
% mask.image = 1/sqrt(size(im,3))*mask.image;

% opA  = load_degmaskoperator('mask',mask,'spectralweights',I_MS_LR.spectralweight,);
opA =load_degmaskoperator('lpfilter',I_MS_LR.KerBlu,'spectralweights',I_MS_LR.spectralweights,'mask',mask,'flag_fft',0);

%% Generation of simulated CASSI acquisition
im=cat(3,I_PAN.data/I_PAN.DynamicRange,I_EXP.data/I_EXP.DynamicRange);
max_value=1; % max_value=I_load{3}.DynamicRange; % If not normalized
y = opA.masking(im);
% y = y+0.001*(255/L)*randn(size(y_ideal));  % Adding very minimal noise 

%% Supporting operators for inversion problem

opL=load_operator('none');
opW=load_operator('CAS8_sym8',size(opA.adjoint(y)));
prox=load_operator('norm_l111');


%% Test for tau Parameter (Most of the times tau=1 is the best choice)

tau_v=[0.1,0.5,1,1.5,2];
rho=1;
lambda=0.01; % Testing value; not necessarily the best one
Nbiter=200;
tol=0; % tol=10^(-6);

x_tau=[]; cost_tau=cell(length(tau_v),1);
for ii=1:length(tau_v)

    fprintf('Evaluation for tau=%.4E\n',tau_v(ii))
    [x_temp,cost_temp]=Loris(y,'opA',opA,'opL',opL,'opW',opW,'prox',prox,'rho',rho,...
        'tol',tol,'Nbiter',Nbiter,'lambda',lambda*max_value,'tau',tau_v(ii),'initx',zeros(size(opA.adjoint(y))));
    
    x_tau=cat(5,x_tau,x_temp);
    cost_tau{ii}=cost_temp{1};
end
x_tau=x_tau*I_MS_LR.DynamicRange;

fprintf('SSIM calculation:\n')
res=zeros(size(x_tau,5),1);

for jj=1:size(x_tau,5)     
    for ii=1:I_GT.size(3)
        res(jj)=res(jj)+ssim(x_tau(:,:,ii,1,jj),I_GT.data(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),max_value);
    end
    res(jj)=res(jj)/size(im,3);
    fprintf('tau: %.4E, SSIM: %6.4f\n',tau_v(jj),res(jj));
end

figure; hold on; title('Primal cost Function'); ylabel('Value'); xlabel('Iteration'); set(gca, 'YScale', 'log')
str_legend=cell(1,size(cost_tau,2));
for jj=1:numel(cost_tau)
    plot(1:numel(cost_tau{jj}),cost_tau{jj}); str_legend{jj}=sprintf('\\tau=%.2f',tau_v(jj));
end
legend(str_legend);

figure; hold on; title('Primal cost Function Variation'); ylabel('Value'); xlabel('Iteration'); set(gca, 'YScale', 'log');
str_legend=cell(1,numel(cost_tau));
for jj=1:numel(cost_tau)
    plot(2:length(cost_tau{jj}),abs(diff(cost_tau{jj},1,1))); str_legend{jj}=sprintf('\\tau=%.2f',tau_v(jj));
end
legend(str_legend);


%% Test for Over-relaxation parameter (Given tau, test for rho=1.5 and 1.9)

lambda=0.001; 
Nbiter=200; tol=0; % Putting tolerance=0 allows a constant number of iterations
tau=1;
rho_v=[1,1.5,1.9,1.99];

x_rho=[]; cost_rho=cell(length(rho_v),1);
for ii=1:length(rho_v)

    [x_temp,cost_temp]=Loris(y,'opA',opA,'opL',opL,'opW',opW,'prox',prox,...
        'rho',rho_v(ii),'tol',tol,'Nbiter',Nbiter,'lambda',lambda*max_value,'tau',tau);
    
    x_rho=cat(5,x_rho,x_temp);
    cost_rho{ii}=cost_temp{1};
end

x_rho=x_rho*I_MS_LR.DynamicRange;

figure; hold on; title('Primal cost Function'); ylabel('Value'); xlabel('Iteration'); set(gca, 'YScale', 'log');
str_legend=cell(1,numel(cost_rho));
for jj=1:numel(cost_rho)
    plot(1:numel(cost_rho{jj}),cost_rho{jj}); str_legend{jj}=sprintf('\\rho=%.2f',rho_v(jj));
end
legend(str_legend);

figure; hold on; title('Primal cost Function Variation'); ylabel('Value'); xlabel('Iteration'); set(gca, 'YScale', 'log')
str_legend=cell(1,size(cost_rho,2));
for jj=1:numel(cost_rho)
    plot(2:numel(cost_rho{jj}),abs(diff(cost_rho{jj},1,1))); str_legend{jj}=sprintf('\\rho=%.2f',rho_v(jj));
end
legend(str_legend);

%% Test for Regularization Parameter (lambda) 
% Note: Algorithm converges very slowly with very low or absent 
% levels of additive Gaussian noise

lambda_v=logspace(-3,-1,3);
tau=1;
rho=1.9;
Nbiter=1000;
tol=0;

[x_lambda,cost_temp]=Loris(y,'opA',opA,'opL',opL,'opW',opW,'prox',prox,...
    'rho',rho,'tau',tau,'lambda',lambda_v*max_value,'init',opA.adjoint(y),'Nbiter',Nbiter,...
    'tol',tol);
x_lambda=x_lambda*I_MS_LR.DynamicRange;


fprintf('SSIM calculation:\n')
res=zeros(size(x_lambda,4),1);

for jj=1:size(x_lambda,4)     
    for ii=1:I_GT.size(3)
        res(jj)=res(jj)+ssim(x_lambda(:,:,ii,jj),I_GT.data(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),max_value);
    end
    res(jj)=res(jj)/size(im,3);
    fprintf('lambda: %5.3f, SSIM: %6.4f\n',lambda_v(jj),res(jj));
end

figure; hold on; title('SSIM Evaluation');
ylabel('Value'); xlabel('Regularization parameter (\lambda)');
plot(lambda_v,res);

[~,idx_lambda]=max(res);

filename=sprintf('%s_%s_WAV',im_tag,mask_tag);
save(fullfile(output_folder,[filename,'.mat']));

figure(1); subplot(241); a=viewimage_outputonly(I_GT.data(:,:,I_GT.Bands_to_display)); imshow(a,[]); title('GT'); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_GT']),'png');
figure(1); subplot(242); a=viewimage_outputonly(I_EXP.data(:,:,I_EXP.Bands_to_display)); imshow(a,[]); title('EXP'); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_EXP']),'png');
figure(1); subplot(243); a=viewimage_outputonly(I_PAN.data(:,:,I_PAN.Bands_to_display)); imshow(a,[]); title('PAN'); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_PAN']),'png');
figure(1); subplot(244); a=viewimage_outputonly(y); imshow(a,[]); title('CASSI Acquisition'); imshow(a,[]); title('PAN'); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_COMP']),'png');
init=opA.adjoint(y);
figure(1); subplot(245); a=viewimage_outputonly(init(:,:,I_MS_LR.Bands_to_display)); imshow(a,[]); title('Initialization'); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_INIT']),'png');
figure(1); subplot(246); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,1)); imshow(a,[]); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(1))); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_INVMIN']),'png');
figure(1); subplot(247); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,idx_lambda)); imshow(a,[]); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(idx_lambda))); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_INVBEST']),'png');
figure(1); subplot(248); a=viewimage_outputonly(x_lambda(:,:,I_MS_LR.Bands_to_display,end)); imshow(a,[]); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(end))); figure(2); imshow(a,[]); saveas(gcf,fullfile(output_folder,[filename,'_INVMAX']),'png');