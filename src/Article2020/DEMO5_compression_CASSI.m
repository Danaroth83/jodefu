clearvars; close all;

% im_tag='Washington_cut256_RGB'; mask='CASSI'; alpha=[];
% im_tag='Janeiro_cut256_RGB'; mask='Bayer'; alpha=[];
% im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; alpha=[];
% im_tag='Janeiro_cut256_all'; mask='period'; alpha=[];
% im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; alpha=[];
% im_tag='Stockholm_cut256_RGB'; mask='Bayer'; alpha=[];
% im_tag='Stockholm_cut256_4'; mask='period'; alpha=[];
% im_tag='Stockholm_cut256_8'; mask='CASSI'; alpha=[];
% im_tag='Janeiro_cut256_8'; mask='CASSI'; alpha=[];
% im_tag='Janeiro_cut256_4'; mask='period'; alpha=[];
% im_tag='Beijing_WV3_WV3_cut256_RGB'; mask='Bayer'; alpha=[];
% im_tag='RdJ_cut256_RGB';
% im_tag='Beijing_WV3_WV3_4';
% im_tag='SanFrancisco_QB_QB_4';
% im_tag='Hobart';

ratio=2;
Nbiter=250; tol=0; % Nbiter=500; tol=10E-06;
alpha=[];
sim_string=0; % If 1, simulated PAN

% preproc='hism'; % preproc='regrnonneg'; % preproc='regr';
% preproc='regravg'; % preproc='regrsum1'; % preproc='none';

%% Non simulated dataset
 
% Tests_image = [1,2,3];
% Tests_kk = [1,2];
% Tests_ii = [1,3];
% Tests_jj = [1,7];

Tests_ii = 0;
Tests_kk = 2;
Tests_jj= 1;
Tests_image = 4;

% Tests_kk=[1,2];
% Tests_ii=1; Tests_jj=[1,7];
% Tests_ii=2; Tests_jj=2;
ra_choice=1.3;

for aa=Tests_image
    if aa==1, im_tag='Beijing_WV3_WV3_4'; end
    if aa==2, im_tag='SanFrancisco_QB_QB_4'; end
    if aa==3, im_tag='Hobart'; end
    if aa==4, im_tag='Washington_cut256_4'; end
    
    for kk=Tests_kk
        if kk==1
            mask='CASSI'; 
            preproc='avg';
            if aa==2, preproc='none'; end
            lambda_v=logspace(-4,-3,5);
        elseif kk==2 
            mask='BinaryTreeU'; 
            preproc='regravg';
            if aa==2, preproc='hism'; end
            lambda_v=logspace(-4,-3,5);
        end

        for ii=Tests_ii
            if ii==0, testtype='nomask'; ra=[]; end
            if ii==1, testtype='default'; ra=[]; end
            if ii==2, testtype='default'; ra=ra_choice; end

            if ii<=2, sim_string=0; SNR=[]; end
            if ii>=3,  sim_string=1; SNR=25; end

            for jj=Tests_jj
                if jj==1, inversion={'TV_c','norm_l221','none'}; end
                if jj==2, inversion={'TV_c','norm_S1l1','none'}; end
                if jj==3, inversion={'TV_u','norm_l221','none'}; end
                if jj==4, inversion={'TV_u','norm_S1l1','none'}; end
                if jj==5, inversion={'TV_s2','norm_l221','none'}; end
                if jj==6, inversion={'TV_s2','norm_S1l1','none'}; end
                if jj==7, inversion={'none','norm_l111','CAS8_sym8'}; end

                fprintf('Mask: %s, Testtype: %s, Sim: %d, Parameters setup: %d, Radius: %.2f\n',mask, testtype,sim_string,jj,ra);

                [I_out,I_acq,mask_out,MR]=wrapper_compressedacquisition('im',im_tag,...
                    'ratio',ratio,'mask',mask,'inv',inversion,'iter',Nbiter,'test',testtype,...
                    'lambda',lambda_v,'preproc',preproc,'tol',tol,...
                    'alpha',alpha,'sim',sim_string,'SNR',SNR,'radius',ra,...
                    'idx_metric','q2n');

            end
        end
    end
end
