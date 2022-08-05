clearvars; close all;

% im_tag='Washington_cut256_4'; mask='period'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Washington_cut256_RGB'; mask='Bayer'; demosaic_list={'WB','ID','IID','ARI2','MLRI2','AP','MSG'};
% im_tag='Janeiro_cut256_RGB'; mask='Bayer'; demosaic_list={'WB','ID','IID','ARI2','MLRI2','AP','MSG'};
% im_tag='Janeiro_cut256_4'; mask='period'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Washington_4'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};

ratio=2;
test_snr = 1;
interpolation='RBF_spline';   % interpolation='WB';
output_folder='compression_classic';

% fusion_list = {'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD','BayesNaive'};
fusion_list = {'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD'};
% fusion_list = {'EXP','PCA','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-HPM','CNMF','BayesNaive'};
% fusion_list = {'EXP', 'GSA', 'BDSD', 'MTF-GLP-HPM', 'MTF-GLP-CBD'};
% fusion_list = {'GSA'};

% testtype = {'default','msonly','nomask'}; % testtype='nodegrad';
testtype = {'default'};

for kk = 1:3

    if kk == 1
        im_tag='Janeiro_RGB';
        mask='Bayer'; 
        demosaic_list={'WB','ID','IID','ARI2','MLRI2','AP','MSG'};
        % demosaic_list = {'ARI2', 'ID'};
    elseif kk==2
        im_tag='Janeiro_4';
        mask='BinaryTreeU'; 
        demosaic_list={'WB','ID','IID','SD','ISD'};
    elseif kk==3
        im_tag='Janeiro_all';
        mask='BinaryTreeU'; 
        demosaic_list={'WB','ID','IID','SD','ISD'};
    end

    for jj=test_snr
    
	    if jj==1, sim_label=0; SNR_db=[]; end
	    if jj==2, sim_label=1; SNR_db=25; end
	    
	    for ii=1:numel(testtype)
	    
	    [I_out,I_acq,mask_out,MR]=wrapper_classic('im',im_tag,...
		    'ratio',ratio,'mask',mask,'interpolation',interpolation,...
		    'fusion',fusion_list,'demosaic',demosaic_list,...
		    'test',testtype{ii},'sim',sim_label,'SNR',SNR_db,...
            'output', output_folder);
		    
	    end
    end
end