clearvars; close all;

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'load_info'),...
        fullfile(project_folder,'scale'),...
        fullfile(project_folder,'visualization'));

%% TAGS FOR MULTIPLATFORM PANSHARPENING    

% im_tag='RdJ_WV3_WV3'; im_tag='RdJ_cut2_HYP_ALI'; im_tag='RdJ_HYP_WV3'; % im_tag='RdJ_ALI_ALI';
% im_tag='Sydney_WV3_WV3'; im_tag='Sydney_HYP_ALI'; % im_tag='Sydney_HYP_WV3';
% im_tag='SaoPaulo_HYP_ALI'; im_tag='SaoPaulo_ALI_ALI'; im_tag='SaoPaulo_IKONOS_IKONOS';  im_tag='SaoPaulo_HYP_IKONOS';

%% TAGS FOR MONOPLATFORM PANSHARPENING
% im_tag='Janeiro'; % im_tag='Janeiro_RGB'; % im_tag='Janeiro_4';
im_tag='Janeiro_cut2';
% im_tag='Tripoli'; % im_tag='Tripoli_RGB'; % im_tag='Tripoli_4';
% im_tag='Washington'; % im_tag='Washington_RGB';
% im_tag='Stockholm'; % im_tag='Stockholm_RGB'; im_tag='Stockholm_4';
% im_tag='Beijing_WV3_WV3'; % im_tag='Beijing_WV3_WV3_RGB';
% im_tag='SanFrancisco_QB_QB'; im_tag='SanFrancisco_cut2_QB_QB_RGB';
% im_tag='Sydney_WV3_WV3';
% im_tag='SaoPaulo_IKONOS_IKONOS'; im_tag='SaoPaulo_cut2_IKONOS_IKONOS';
% im_tag='Vancouver'; im_tag='Moffett'; % im_tag='Moffett_cut2'; im_tag='Moffett_ovlp';  im_tag='Moffett_all';
% im_tag='Pavia';
% im_tag='Rome'; % im_tag='China'; % im_tag='Pleiades'; % im_tag='Pleiades_cut2';
% im_tag='Toulouse'; im_tag='Toulouse_cut2';
% im_tag='Indianapolis'; % im_tag='Indianapolis_cut2';
% im_tag='Hobart'; im_tag='Hobart_cut2';
% im_tag='Sofia'; % im_tag='Sofia_ALI_ALI'; % im_tag='Sofia_cut2_ALI_ALI_VNIR';
% im_tag='Sudbury'; % im_tag='Sudbury_ALI_ALI'; % im_tag='Sofia_cut2_ALI_ALI_VNIR';
% im_tag='Fuji'; % im_tag='Fuji_ALI_ALI'; % im_tag='Fuji_ALI_ALI_VNIR';
% im_tag='Paris'; im_tag='Chris';

I_load=Load_Dataset_Pansharpening(im_tag,'test','FR','request',{'MS_LR','PAN','REF'});

for ii=1:numel(I_load)
    I=I_load{ii};
    fprintf('type:   %s\n',I.type);
    fprintf('label:  %s\n',I.label);
    fprintf('sensor: %s\n',I.sensor);
    fprintf('GSD:    %.3f\n',I.GSD);
    fprintf('sizes:  %d %d %d\n',size(I.data,1),size(I.data,2),size(I.data,3));
    fprintf('ratio:  %d\n',I.ratio);
    viewimage(I.data(:,:,I.Bands_to_display));
    fprintf('\n');
end

out=dftregistration(fft2(I_load{2}.data),fft2(I_load{3}.data(:,:,1)),11);
fprintf('Disalignment: %f,%f\n',out(3),out(4));

a=show_RGB(I_load{3});