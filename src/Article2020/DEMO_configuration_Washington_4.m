clearvars; close all;

% im_tag='Washington_cut256_RGB'; mask='Bayer'; alpha=[];
% im_tag='Janeiro_cut256_RGB'; mask='Bayer'; alpha=[];
% im_tag='Janeiro_cut256_all'; mask='period'; alpha=[];
% im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; alpha=[];
% im_tag='Janeiro_cut256_4'; mask='period'; alpha=[];
% im_tag='Beijing_WV3_WV3_cut256_RGB'; mask='Bayer'; alpha=[];
im_tag='Washington_cut256_4'; mask='BinaryTreeU'; alpha=[];
ratio=2;
Nbiter=250; tol=0; % Nbiter=500; tol=10E-06;


% preproc='hism';
% preproc='regrnonneg';
% preproc='regr';
% preproc='regravg';
% preproc='regrsum1';
% preproc='none'; % lambda_v=logspace(-2.5,-1.5,5);

%% Non simulated dataset

% Tests_ii=1:8;
% Tests_ii=[1,3,5,6,8];
% Tests_ii=[2,4,7];
% Tests_ii=[3,4,6,7];
Tests_ii=[2,4,5];
% Tests_ii=1:8;
%Tests_ii=[2,4];


% Tests_jj=1:6;
Tests_jj=4;
% Tests_jj=1;

ra_choice=1.5;
ra_choice2=1.3;

for ii=Tests_ii   
    if ii==1, testtype='default';  lambda_v=logspace(-3.2,-2.2,5); ra=[]; end
    if ii==2, testtype='default';  lambda_v=logspace(-3.2,-2.2,5); ra=ra_choice; end
    if ii==3, testtype='msonly'; lambda_v=logspace(-3.2,-2.2,5); ra=[]; end
    if ii==4, testtype='msonly'; lambda_v=logspace(-3.2,-2.2,5); ra=ra_choice2; end    
    if ii==5, testtype='nomask'; lambda_v=logspace(-4,-3,5); ra=[]; end
    
    if ii==6, testtype='default';  lambda_v=logspace(-3.5,-2.5,5); ra=[]; end
    if ii==7, testtype='default';  lambda_v=logspace(-3.5,-2.5,5); ra=ra_choice; end
    if ii==8, testtype='nomask'; lambda_v=logspace(-4,-3,5); ra=[]; end
    
    %lambda_v=logspace(-3.5,-2.5,5);
    
    if ii<=5, sim_string=0; SNR=[]; preproc='hism'; end
    if ii>5,  sim_string=1; SNR=25; preproc='hism'; end
        
    for jj=Tests_jj
        if jj==1, inversion={'TV_c','norm_l221','none'}; end
        if jj==2, inversion={'TV_c','norm_S1l1','none'}; end
        if jj==3, inversion={'TV_u','norm_l221','none'}; end
        if jj==4, inversion={'TV_u','norm_S1l1','none'}; end
        if jj==5, inversion={'TV_s2','norm_l221','none'}; end
        if jj==6, inversion={'TV_s2','norm_S1l1','none'}; end

        fprintf('Testtype: %s, Sim: %d, Parameters setup: %d, Radius: %.2f\n',testtype,sim_string,jj,ra);

        [I_out,I_acq,mask_out,MR]=wrapper_compressedacquisition('im',im_tag,...
            'ratio',ratio,'mask',mask,'inv',inversion,'iter',Nbiter,'test',testtype,...
            'lambda',lambda_v,'preproc',preproc,'tol',tol,'alpha',alpha,'sim',sim_string,...
            'SNR',SNR,'radius',ra);
        
    end
end
