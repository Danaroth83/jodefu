clearvars; close all;

current_folder=fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'..','jodefu'));

output_folder = 'formation_classic';

% im_tag='Washington_cut256_4'; mask='period'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Washington_cut256_RGB'; mask='Bayer'; demosaic_list={'WB','ID','IID','ARI2','MLRI2','AP','MSG'};
% im_tag='Janeiro_cut256_RGB'; mask='Bayer'; demosaic_list={'WB','ID','IID','ARI2','MLRI2','AP','MSG'};
% im_tag='Janeiro_cut256_4'; mask='period'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Stockholm_cut256_all'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};
im_tag='Washington_4'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};
% im_tag='Janeiro_cut256_all'; mask='BinaryTreeU'; demosaic_list={'WB','ID','IID','SD','ISD'};

ratio=2;

interpolation='RBF_spline';    % interpolation='WB';

% fusion_list={'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD','BayesNaive'};
fusion_list={'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD'};
%fusion_list={'EXP','PCA','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-HPM','CNMF','BayesNaive'};

testtype={'default','msonly','nomask'}; % testtype='nodegrad';
tests_sim = 1; % 1= Non simulated PAN, 2 = PAN simulated from GT (SNR = 25 dB)

for jj=tests_sim
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