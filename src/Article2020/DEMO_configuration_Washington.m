clearvars; close all;

im_tag='Washington_cut256_4'; mask='period'; alpha=[];
% im_tag='Washington_cut256_RGB'; mask='Bayer'; alpha=[];
ratio=2;
% inversion='TV';
Nbiter=250; tol=0; % Nbiter=500; tol=10E-06;

% inversion={'TV_c','norm_l221','none'};
% inversion={'TV_c','norm_S1l1','none'};
% inversion={'TV_s2','norm_l221','none'};
% inversion={'TV_s2','norm_S1l1','none'};
% inversion={'TV_u','norm_l221','none'};
inversion={'TV_u','norm_S1l1','none'};

% testtype='default';
% testtype='msonly'; % Note: This setup is completely deterministic as no noise is added. This means that low lambda should work the best as the function is fully invertible
% testtype='nodegrad';
% testtype='nomask';

% preproc='hism';
% preproc='regrnonneg';
% preproc='regr';
% preproc='regravg';
% preproc='regrsum1';
% preproc='none'; % lambda_v=logspace(-2.5,-1.5,5);

%% Non simulated dataset

%sim=0; SNR=[]; preproc='hism';
%testtype='default';  lambda_v=logspace(-3,-2,5);
% testtype='msonly';  lambda_v=logspace(-3.5,-2.5,5); % Note: This setup is completely deterministic as no noise is added. This means that low lambda should work the best as the function is fully invertible
% testtype='nodegrad'; lambda_v=logspace(-3,-2,5);
% testtype='nomask'; lambda_v=logspace(-4,-3,5);

%% Simulated dataset

sim=1; SNR=25; preproc='hism';
testtype='default';  lambda_v=logspace(-3,-2,5);
% testtype='msonly';  lambda_v=logspace(-3.5,-2.5,5); % Note: This setup is completely deterministic as no noise is added. This means that low lambda should work the best as the function is fully invertible
% testtype='nodegrad'; lambda_v=logspace(-3,-2,5);
% testtype='nomask'; lambda_v=logspace(-3,-2,5);


% im_tag='Washington_cut256_4'; mask='random'; alpha=0.75;
% preproc='hism'; lambda_v=logspace(-2,-1,5);

% im_tag='Washington_cut256_4'; mask='random'; alpha=0.5;
% preproc='hism'; lambda_v=logspace(-2,-1,5);

[I_out,I_acq,mask,MR]=wrapper_compressedacquisition('im',im_tag,...
    'ratio',ratio,'mask',mask,'inv',inversion,'iter',Nbiter,'test',testtype,...
    'lambda',lambda_v,'preproc',preproc,'tol',tol,'alpha',alpha,'sim',sim,...
    'SNR',SNR);

