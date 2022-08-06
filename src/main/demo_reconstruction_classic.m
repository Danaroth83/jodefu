clearvars; close all;

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'jodefu'));

ratio=2; % Scale ratio
interpolation = 'RBF_spline'; % PAN interpolation
fusion_list={'EXP','GSA','BDSD','ATWT','MTF-GLP-HPM','MTF-GLP-CBD'}; % Fusion methods

%% Image reconstruction; JoDeFu results
output_folder = 'reconstruction_classic'; % Output folder
image_list = 1:3; % Janeiro with 1 = 3 bands, 2 = 4 bands, 3 = 8 bands
formation_list = 1; % 0 = Spa./spe. degradation, 1 = MRCA, 2 = Mosaicing 

for kk=image_list
    preproc='regravg';
    if kk==1, im_tag = 'Janeiro_3'; mask = 'Bayer'; end
    if kk==2, im_tag = 'Janeiro_4'; mask = 'BinaryTreeU'; end
    if kk==3, im_tag = 'Janeiro_8'; mask = 'BinaryTreeU'; end
    
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