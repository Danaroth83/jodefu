function demo_formation_classic(varargin)

if numel(varargin) >= 1, im_tag = varargin{1}; else, im_tag = 'Washington_4'; end

%clearvars; 
close all;
fprintf('==== Image formation tests (Classic methods) ====\n');

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'jodefu'));

%% Image formation - Classic models

output_folder = 'formation_classic'; % Output folder
ratio=2; % Scale ratio
% im_tag = 'Washington_4'; % Image tag
mask = 'BinaryTreeU'; % Mask label
interpolation='RBF_spline';    % PAN interpolation
demosaic_list = {'WB', 'ID', 'IID', 'SD', 'ISD'}; % Demosaic method list
fusion_list={'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD'}; % Fusion list
formation_list = 1:3; % 1=Spa./spe. degradation, 2=MRCA, 3=Mosaicing

for jj=formation_list
    sim_label=0; SNR_db=[];
    if jj==1, testtype = 'nomask'; end
    if jj==2, testtype = 'default'; end
    if jj==3, testtype = 'msonly'; end
    if jj==4, testtype = 'default'; sim_label=1; SNR_db=25; end
	
	[I_out,I_acq,mask_out,MR]=wrapper_classic('im',im_tag,...
		'ratio',ratio,'mask',mask,'interpolation',interpolation,...
		'fusion',fusion_list,'demosaic',demosaic_list,...
		'test',testtype,'sim',sim_label,'SNR',SNR_db,...
        'output', output_folder);
end