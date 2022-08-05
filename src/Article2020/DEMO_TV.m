%% DENOISING WITH TOTAL VARIATION
% 
% Author:
% Daniele Picone
% v1.0
% 
% Description:
% Inversion of a noisy acquisition via anisotropic Total Variation
% The algorithm solves the problem of denoising by minimizing the cost
% function
% f(x)= (1/2) ||x-y||_2^2 + lambda ||Lx||_2,1
% where:
%  - y is the noisy image
%  - L is the gradient in the horizontal and vertical direction
%  - ||.||_2 is a norm l_2 over all dimensions
%  - ||.||_2,1 is a norm l_2 over the gradients and l_1 over the other
%       dimensions


%% Adding support path

clearvars; close all;
support_folder=fullfile('..','..','support');
addpath(fullfile(support_folder,'Operator'),...
        fullfile(support_folder,'Quality_indices'),...
        fullfile(support_folder,'Inversion_Solver'),...
        fullfile(support_folder,'Visualization'));

%% Generation of noisy image
L=255; % If =1, no normalization
im=double(imread('peppers.png'))/L;
y=im+0.1*(255/L)*randn(size(im));

lambda=0.05:0.01:0.15;

opA=load_operator('none');
opg=load_operator('Norm_l211');
opL=load_operator('TVc');

[x,cost]=Loris(y,'Nbiter',250,'opA',opA,'opL',opL,'reg',opg,...
               'lambda',lambda*(255/L),'rho',1.9,'tau',1,'tol',10E-6,...
               'init',opA.adjoint(y));

res=zeros(1,length(lambda));
for jj=1:size(x,4)     
    for ii=1:size(im,3)
        res(jj)=res(jj)+ssim(x(:,:,ii,jj),im(:,:,ii),[0.01,0.03],fspecial('gaussian', 11, 1.5),255/L);
    end
    res(jj)=res(jj)/size(im,3);
    fprintf('lambda: %.4E, SSIM: %6.4f\n',lambda(jj)*255/L,res(jj));
end

figure; hold on; title('SSIM Evaluation');
ylabel('Value'); xlabel('Regularization parameter (\lambda)');
plot(lambda,res);

[~,idx_lambda]=max(res);

subplot(221); imshow(im,[]); title('Original');
subplot(222); imshow(y,[]); title('Noisy Image');
subplot(223); imshow(opA.adjoint(y),[]); title('Inversion (First \lambda)');
subplot(224); imshow(x(:,:,:,idx_lambda),[]); title('Inversion (Best)');