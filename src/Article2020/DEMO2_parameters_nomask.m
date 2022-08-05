clearvars; close all;

output_main=fullfile('..','..','data','output');
output_fol='test_parameters';
ratio=2;
Nbiter=250; tol=0; % Nbiter=500; tol=10E-06;
% Nbiter=250; tol=0;
sim_string=0; %if 1, simulated dataset

% test={'lambda','norm','blur','tv','preproc'};
test={'best'};
% test={'preproc'};

% test={'lambda'};
Tests_kk=1;

idx_best='Q2^n'; %idx_best='SSIM'; %idx_best='PSNR'; idx_best='ERGAS';

% Tests_ii=4;
% if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; ra_choice=1.5; lambda_v=1.00E-03; end
% if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; ra_choice=1.5; lambda_v=1.00E-03; end
% if kk==3, im_tag='Washington_cut256_4'; mask='period'; ra_choice=1.5; lambda_v=1.00E-03; end
% if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; ra_choice=1.5; lambda_v=1.78E-03; end
% if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; ra_choice=1.5; lambda_v=2.00E-03; end
% if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; ra_choice=1.5; lambda_v=3.16E-03; end

for tt=1:numel(test)

    for kk=Tests_kk
        
        ra=1; tv={'TV_u'}; normtype={'norm_S1l1'}; preproc={'regravg'};
        testtype={'nomask'};
        % testtype={'mixed'};
        demosaic_list={'ARI2','MLRI2'};

        if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; lambda_v=1.78E-03; end %1.93E-03
        if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; lambda_v=1.78E-03; end %1.93E-03
        if kk==3, im_tag='Washington_cut256_4'; mask='period'; lambda_v=5.62E-04; end %5.62E-04;
        if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; lambda_v=1.00E-03; end %7.20E-04;
        if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; lambda_v=1.12E-03; end %9.21E-04
        if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; lambda_v=1.78E-03; end 
        if kk==7, im_tag='Beijing_cut256_WV3_WV3_RGB'; mask='Bayer'; lambda_v=5.62E-04; end
        if kk==8, im_tag='Beijing_cut256_WV3_WV3_4'; mask='period'; lambda_v=5.62E-04; end
        
        if strcmpi(test{tt},'best')
            lambda_v=logspace(-2.5,-2,5); 
            var=lambda_v; 
            x_label='Regularization Parameter (\lambda)'; 
            output_fol=fullfile('test_parameters','best');
        else
            output_fol=fullfile('test_parameters');
            if strcmpi(test{tt},'lambda'), lambda_v=logspace(-4,-3,20); var=lambda_v; x_label='Regularization Parameter (\lambda)'; output_fol=fullfile('test_parameters','lambda'); end
            if strcmpi(test{tt},'blur'), ra=[1,1.3:0.1:2]; var=ra; x_label='Blur Diameter'; else, ra=1; end
            if strcmpi(test{tt},'tv'), tv={'TV_c','TV_u','TV_s2','TV_s3'}; var=tv; x_label='Total Variation'; else, tv={'TV_c'}; end
            if strcmpi(test{tt},'norm'), normtype={'norm_l221','norm_S1l1','norm_l211','norm_linf11','norm_l111','norm_Sinfl1'}; var=normtype; x_label='norm'; else, normtype={'norm_l221'}; end
            if strcmpi(test{tt},'preproc'), preproc={'none','hism','regrnonneg','regr','regravg','regrsum1'}; var=preproc; x_label='Pre-Processing'; else, preproc={'hism'}; end
            if strcmpi(test{tt},'methods'), testtype={'default','msonly','nomask'}; var=testtype; x_label='Test Type'; else, testtype={'default'}; end
        end
            
        %% Non simulated dataset
        
        fprintf('Dataset: %s\n',im_tag);
        MR_out=[];
        for ii=numel(testtype)   
            alpha=[];
            if sim_string==0, SNR=[]; else, SNR=25; end 
            for jjjj=1:numel(preproc)
                for ll=1:numel(ra)
                    for jj=1:numel(normtype)
                        for jjj=1:numel(tv)      

                            inversion={tv{jjj},normtype{jj},'none'};
                            fprintf('Testtype: %s, Sim: %d, Norm: %s, TV: %s, Radius: %.2f, Preproc: %s\n',testtype{ii},sim_string,normtype{jj},tv{jjj},ra(ll),preproc{jjjj});

                            [I_out,I_acq,mask_out,MR]=wrapper_compressedacquisition(...
                                'im',im_tag,'ratio',ratio,'mask',mask,'inv',inversion,...
                                'iter',Nbiter,'test',testtype{ii},'lambda',lambda_v,...
                                'preproc',preproc{jjjj},'tol',tol,'alpha',alpha,'sim',sim_string,...
                                'SNR',SNR,'radius',ra(ll),'idx_metric','q2n',...
                                'output',output_fol,'demosaic',demosaic_list);

                            MR_out=cat(2,MR_out,MR.data);
                            
                        end       
                    end
                end
                %% Plot Parameter Graph
                figure; 
                if iscell(var)
                    xticks(1:numel(var)); xticklabels(var); var_out=1:numel(var);
                else
                    var_out=var;
                end
                for iii=1:numel(MR.qindex), if strcmpi(MR.qindex{iii},idx_best), idx_best=iii; end, end
                plot(var_out,MR_out(idx_best,:)); xlabel(x_label); ylabel(MR.qindex{idx_best});
                str_iter=sprintf('i%d',Nbiter);
                filename=fullfile(output_main,output_fol,[im_tag,'_',test{tt},'_',testtype{ii},sprintf('_i%d_',Nbiter),preproc{jjjj}]);
                drawnow;
                saveas(gcf,[filename,'.fig']);
                saveas(gcf,[filename,'.png']);
            end
        end
    end
end