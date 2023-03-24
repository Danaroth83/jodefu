function demo_reconstruction_classic(varargin)

if numel(varargin) >= 1, im_tag = varargin{1}; else, im_tag = 'Janeiro'; end

% clearvars; 
close all;
fprintf('== Image reconstruction tests (Classic methods) ==\n');

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'jodefu'));

%% Image reconstruction - Classic algorithms
output_folder = 'reconstruction_classic'; % Output folder
if strcmpi(im_tag, 'Janeiro')
    image_list = 1:3; % Janeiro with 1 = 3 bands, 2 = 4 bands, 3 = 8 bands
elseif strcmpi(im_tag, 'Stockholm')
    image_list = 4:6; % Stockholm with 1 = 3 bands, 2 = 4 bands, 3 = 8 bands
else
    error('Available datasets: [Janeiro, Stockholm]');
end
formation_list = 1; % 0 = Spa./spe. degradation, 1 = MRCA, 2 = Mosaicing 

ratio=2; % Scale ratio
interpolation = 'RBF_spline'; % PAN interpolation
fusion_list={'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD'}; % Fusion methods

for kk=image_list
    preproc='regravg';
    if kk==1, im_tag = 'Janeiro_3'; mask = 'Bayer'; end
    if kk==2, im_tag = 'Janeiro_4'; mask = 'BinaryTreeU'; end
    if kk==3, im_tag = 'Janeiro_8'; mask = 'BinaryTreeU'; end
    if kk==4, im_tag = 'Janeiro_3'; mask = 'Bayer'; end
    if kk==5, im_tag = 'Janeiro_4'; mask = 'BinaryTreeU'; end
    if kk==6, im_tag = 'Janeiro_8'; mask = 'BinaryTreeU'; end
    
    if kk==1
        demosaic_list = {'WB', 'ID', 'IID', 'ARI2', 'MLRI2', 'AP', 'MSG'};
    else
        demosaic_list = {'WB', 'ID', 'IID', 'SD', 'ISD'};
    end
        
    for ii=formation_list
        sim_string=0; SNR = [];
        if ii==0, testtype='nomask'; end
        if ii==1, testtype='default'; end
        if ii==2, testtype='msonly'; end
        if ii==3, testtype='default'; sim_string=1; SNR=25; end

        [I_out,I_acq,mask_out,MR]=wrapper_classic('im', im_tag,...
            'ratio', ratio, 'mask', mask, 'interpolation', interpolation,...
            'fusion', fusion_list, 'demosaic', demosaic_list,...
            'test', testtype, 'sim', sim_string, 'SNR', SNR,...
            'output', output_folder);

    end
end