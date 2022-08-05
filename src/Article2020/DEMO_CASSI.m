%% TEST FOR CASSI SOLVER
%
% Author:
% Daniele Picone
% v1.0
%
% Description:
% This script allows to solve the inversion problem of a simulated
% acquisition with a CASSI device, according to the procedure described in
% [1]. The testing procedure follows the suggestions of [2]
% 
% References:
% [1] Arce G. R., Brady D. J., Carin L., Arguello H. and Kittle D. S., 
% "Compressive Coded Aperture Spectral Imaging: An Introduction", IEEE SPM,
%  2014, 31, 105-115
% [2] Condat L.,  Kitahara D., Contreras A. and Hirabayashi A., "Proximal
% splitting algorithms: Relax them all!" (preprint)



%% Adding support path

clearvars; close all;
rng('default');  % For reproductible results
current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'operator'));
addpath(fullfile(project_folder,'quality_indices'));
addpath(fullfile(project_folder,'inversion_solver'));
addpath(fullfile(project_folder,'mosaic'));


%% Loading Image and Generating Acquisition

L=255; % 1= no normalization; 255= normalization
im=double(imread('peppers.png'))/L;
im=im(size(im,1)/2+(-127:128),size(im,2)/2+(-127:128),:);

%% Setting Operator for masking operation
mask = load_mask('mask','CASSI','sizes',size(im));
% mask.image = 1/sqrt(size(im,3))*mask.image;
opA  = load_degmaskoperator('mask',mask);

%% Generation of simulated CASSI acquisition
y = opA.direct(im);
% y = y+0.001*(255/L)*randn(size(y_ideal));  % Adding very minimal noise 

%% Supporting perators for inversion problem

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
        'tol',tol,'Nbiter',Nbiter,'lambda',lambda*(255/L),'tau',tau_v(ii),'initx',zeros(size(opA.adjoint(y))));
    
    x_tau=cat(5,x_tau,x_temp);
    cost_tau{ii}=cost_temp{1};
end


fprintf('SSIM calculation:\n')
res=zeros(size(x_tau,5),1);

for jj=1:size(x_tau,5)     
    for ii=1:size(im,3)
        res(jj)=res(jj)+ssim(x_tau(:,:,ii,1,jj),im(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),255/L);
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

lambda=0.01; 
Nbiter=200; tol=0; % Putting tolerance=0 allows a constant number of iterations
tau=1;
rho_v=[1,1.5,1.9,1.99];

x_rho=[]; cost_rho=cell(length(rho_v),1);
for ii=1:length(rho_v)

    [x_temp,cost_temp]=Loris(y,'opA',opA,'opL',opL,'opW',opW,'prox',prox,...
        'rho',rho_v(ii),'tol',tol,'Nbiter',Nbiter,'lambda',lambda*(255/L),'tau',tau);
    
    x_rho=cat(5,x_rho,x_temp);
    cost_rho{ii}=cost_temp{1};
end

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

lambda_v=logspace(-4,-2,3);
tau=2;
rho=1;
Nbiter=2000;
tol=0;

[x_lambda,cost_temp]=Loris(y,'opA',opA,'opL',opL,'opW',opW,'prox',prox,...
    'rho',rho,'tau',tau,'lambda',lambda_v*(255/L),'init',opA.adjoint(y),'Nbiter',Nbiter,...
    'tol',tol);

fprintf('SSIM calculation:\n')
res=zeros(size(x_lambda,4),1);

for jj=1:size(x_lambda,4)     
    for ii=1:size(im,3)
        res(jj)=res(jj)+ssim(x_lambda(:,:,ii,jj),im(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),255/L);
    end
    res(jj)=res(jj)/size(im,3);
    fprintf('lambda: %5.3f, SSIM: %6.4f\n',lambda_v(jj),res(jj));
end


figure; hold on; title('SSIM Evaluation');
ylabel('Value'); xlabel('Regularization parameter (\lambda)');
plot(lambda_v,res);

[~,idx_lambda]=max(res);

figure;
subplot(231); imshow(im,[]); title('Original');
subplot(232); imshow(y,[]); title('CASSI Acquisition');
subplot(233); imshow(opA.adjoint(y),[]); title('Initialization');
subplot(234); imshow(x_lambda(:,:,:,1),[]); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(1)));
subplot(235); imshow(x_lambda(:,:,:,idx_lambda),[]); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(idx_lambda)));
subplot(236); imshow(x_lambda(:,:,:,end),[]); title(sprintf('Inversion (\\lambda=%.2E)',lambda_v(end)));
