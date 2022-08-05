clearvars; close all;

ratio=2;
Nbiter=250; tol=0; % Nbiter=500; tol=10E-06;
% Nbiter=250; tol=0;
sim_string=0; %if 1, simulated dataset


% test={'lambda','norm','blur','tv','preproc'};
%test={'best'};
% test={'lambda2'};
% test={'preproc'};

% test={'norm','blur','tv','preproc','lambda'};
% test={'blur','tv','norm','lambda'};
% test={'norm'};
% test={'best2'};
% test={'blur2'};
% test={'lambda'};
% test={'base', 'lambda', 'tv', 'blur', 'blur2', 'norm', 'lambda2'};
test = {'best'};
Tests_kk=9;
idx_best='SSIM'; %idx_best='Q2^n'; %idx_best='PSNR'; idx_best='ERGAS';
output_fol='test_parameters_new';


% Tests_ii=4;
% if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; ra_choice=1.5; lambda_v=1.00E-03; end
% if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; ra_choice=1.5; lambda_v=1.00E-03; end
% if kk==3, im_tag='Washington_cut256_4'; mask='period'; ra_choice=1.5; lambda_v=1.00E-03; end
% if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; ra_choice=1.5; lambda_v=1.78E-03; end
% if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; ra_choice=1.5; lambda_v=2.00E-03; end
% if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; ra_choice=1.5; lambda_v=3.16E-03; end


for tt=1:numel(test)

    for kk=Tests_kk
        
        ra=1.3; tv={'TV_c'}; normtype={'norm_l221'}; testtype={'default'}; preproc={'regravg'};

        if kk==1, im_tag='Washington_cut256_RGB'; mask='Bayer'; lambda_v=1.93E-03; tv={'TV_s2'}; end %lambda_v=1.78E-03;
        if kk==2, im_tag='Janeiro_cut256_RGB'; mask='Bayer'; lambda_v=1.93E-03; end %lambda_v=1.78E-03;
        if kk==3, im_tag='Washington_cut256_4'; mask='period'; lambda_v=5.62E-04; end %1.00E-03;
        if kk==4, im_tag='Janeiro_cut256_4'; mask='period'; lambda_v=7.20E-04; end %7.20E-04;
        if kk==5, im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; lambda_v=9.21E-04; end %1.12E-03;
        if kk==6, im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; lambda_v=9.21E-04; end %1.78E-03;
        if kk==7, im_tag='Beijing_cut256_WV3_WV3_RGB'; mask='Bayer'; lambda_v=1.92E-03; tv={'TV_c'}; ra=1.5; end %5.62E-04; 
        if kk==8, im_tag='Beijing_cut256_WV3_WV3_4'; mask='period'; lambda_v=9.21E-04; tv={'TV_c'}; ra=1.5; end %5.62E-04;
        if kk==9, im_tag='Beijing_WV3_WV3_4'; mask='BinaryTreeU'; lambda_v=1E-03; tv={'TV_c'}; ra=1.5; end %5.62E-04;
        
        output_fol = fullfile(output_fol, im_tag);

        if strcmpi(test{tt},'best')
            lambda_v=0.001;
            tv = {'TV_s2'};
            ra = 1.3;
            normtype = {'norm_S1l1'};
            var=lambda_v; 
            x_label='Regularization Parameter (\lambda)'; 
            output_fol=fullfile(output_fol,'best');
        elseif strcmpi(test{tt},'best2')
            % lambda_v=logspace(-4,-2.5,15); 
            var=lambda_v; 
            x_label='Regularization Parameter (\lambda)'; 
            output_fol=fullfile(output_fol,'best2');
        else
            output_fol=fullfile(output_fol, im_tag, test{tt});
            if strcmpi(test{tt},'base'), var=lambda_v; x_label='Regularization Parameter (\lambda)'; output_fol=fullfile(output_fol,'base'); end
            if strcmpi(test{tt},'lambda2'), lambda_v=logspace(-6,-1,2); var=lambda_v; x_label='Regularization Parameter (\lambda)';
            elseif strcmpi(test{tt},'lambda'), lambda_v=logspace(-4,-2.5,15); var=lambda_v; x_label='Regularization Parameter (\lambda)'; else, lambda_v=1E-03;  end
            if strcmpi(test{tt},'blur2'), ra=3; var=ra; x_label='Blur Diameter'; output_fol=fullfile('test_parameters','extra');
            elseif strcmpi(test{tt},'blur'), ra=[1,1.3:0.1:2]; var=ra; x_label='Blur Diameter'; else, ra=1; end
            if strcmpi(test{tt},'tv'), tv={'TV_c','TV_u','TV_s2','TV_s3'}; var=tv; x_label='Total Variation'; else, tv={'TV_c'}; end
            if strcmpi(test{tt},'norm'), normtype={'norm_l221','norm_S1l1','norm_l211','norm_linf11','norm_l111','norm_Sinfl1'}; var=normtype; x_label='norm'; else, normtype={'norm_l221'}; end
            if strcmpi(test{tt},'preproc'), preproc={'none','hism','regrnonneg','regr','regravg','regrsum1'}; var=preproc; x_label='Pre-Processing'; else, preproc={'regravg'}; end
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
                                'output',output_fol);

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