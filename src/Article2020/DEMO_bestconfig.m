clearvars; close all;

output_fol='test_moreiter';
ratio=2;
Nbiter=250; tol=0; % Nbiter=500; tol=10E-06;

test_lambda=0; % If 1, tests on all valeus of lambda
%Tests_kk=1:6; %Choice Dataset
Tests_kk=3;
% testtype='default'; 
testtype='msonly'; 
% testtype='nomask';
Tests_jj=4;  % Inversion parameters


% preproc='hism';
% preproc='regrnonneg';
% preproc='regr';
% preproc='regravg';
% preproc='regrsum1';
% preproc='none';

for kk=Tests_kk
    if strcmpi(testtype,'nomask')
        Tests_ii=5; ra_choice=[];
        if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; lambda_v=5.62E-04; end
        if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; lambda_v=1.00E-03; end
        if kk==3, im_tag='Washington_cut256_4'; mask='period'; lambda_v=1.78E-04; end
        if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; lambda_v=3.16E-04; end
        if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; lambda_v=5.62E-04; end
        if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; lambda_v=3.16E-04; end
    elseif strcmpi(testtype,'msonly')
        Tests_ii=4; ra_choice=1.5;
        if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; lambda_v=1.00E-03; end
        if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; lambda_v=1.00E-03; end
        if kk==3, im_tag='Washington_cut256_4'; mask='period'; lambda_v=1.00E-03; end
        if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; lambda_v=1.78E-03; end
        if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; lambda_v=2.00E-03; end
        if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; lambda_v=3.16E-03; end
    else
        Tests_ii=2; ra_choice=1.3;
        if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; lambda_v=1.78E-03; end
        if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; lambda_v=1.78E-03; end
        if kk==3, im_tag='Washington_cut256_4'; mask='period'; lambda_v=5.62E-04; end
        if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; lambda_v=1.00E-03; end
        if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; lambda_v=1.12E-03; end
        if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; lambda_v=1.78E-03; end
    end
    %% Non simulated dataset
    if test_lambda==1, lambda_v=logspace(-4,-2,20); end

    fprintf('Dataset: %s\n',im_tag);
    for ii=Tests_ii   
        if ii==1, testtype='default'; ra=[]; end
        if ii==2, testtype='default'; ra=ra_choice; end
        if ii==3, testtype='msonly'; ra=[]; end
        if ii==4, testtype='msonly'; ra=ra_choice; end    
        if ii==5, testtype='nomask'; ra=[]; end

        if ii==6, testtype='default'; ra=[]; end
        if ii==7, testtype='default'; ra=ra_choice; end
        if ii==8, testtype='nomask'; ra=[]; end

        if ii<=5, sim_string=0; SNR=[]; preproc='hism'; alpha=[]; end
        if ii>5,  sim_string=1; SNR=25; preproc='hism'; alpha=[]; end

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
                'SNR',SNR,'radius',ra,'output',output_fol);

        end
    end
end
