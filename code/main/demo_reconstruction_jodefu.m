function demo_reconstruction_jodefu(varargin)

if numel(varargin) >= 1, im_tag = varargin{1}; else, im_tag = 'Janeiro'; end

% clearvars; 
close all;
fprintf('====== Image reconstruction tests (JoDeFu) ======\n');

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'jodefu'));

%% Image reconstruction; JoDeFu algorithm
output_folder = 'reconstruction_jodefu'; % Output folder
if strcmpi(im_tag, 'Janeiro')
    image_list = 1:3; % Janeiro with 1 = 3 bands, 2 = 4 bands, 3 = 8 bands
elseif strcmpi(im_tag, 'Stockholm')
    image_list = 4:6; % Stockholm with 1 = 3 bands, 2 = 4 bands, 3 = 8 bands
else
    error('Available datasets: [Janeiro, Stockholm]');
end
formation_list = 1; % 0 = Spa./spe. degradation, 1 = MRCA, 2 = Mosaicing
reconstruction_list = [0, 4, 8]; % 0 = jodefu v1, 4&8 = jodefu v2 with UTV
lambda_v = 0.001:0.001:0.002; % Normalized regularization parameter

ratio=2; % Scale ratio
Nbiter=250; % Number of algorithm iterations
tol=0; % cost function tolerance (if 0, stop when iteration = Nbiter)
d_b_choice=1.4; % Blurring diameter for jodefu v2

for kk=image_list
    preproc='regravg';
    if kk==1, im_tag = 'Janeiro_3'; mask = 'Bayer'; end
    if kk==2, im_tag = 'Janeiro_4'; mask = 'BinaryTreeU'; end
    if kk==3, im_tag = 'Janeiro_8'; mask = 'BinaryTreeU'; end
    if kk==4, im_tag = 'Stockholm_3'; mask = 'Bayer'; end
    if kk==5, im_tag = 'Stockholm_4'; mask = 'BinaryTreeU'; end
    if kk==6, im_tag = 'Stockholm_8'; mask = 'BinaryTreeU'; end

    for ii=formation_list
        sim_string=0; SNR = [];
        if ii==0, testtype='nomask'; end
        if ii==1, testtype='default'; end
        if ii==2, testtype='msonly'; end
        if ii==3, testtype='default'; sim_string=1; SNR=25; end

        for jj=reconstruction_list
            if jj==0, inversion={'TV_c','norm_l221','none'}; d_b = []; end
            if jj==1, inversion={'TV_c','norm_l221','none'}; d_b = d_b_choice; end
            if jj==2, inversion={'TV_c','norm_S1l1','none'}; d_b = d_b_choice; end
            if jj==3, inversion={'TV_u','norm_l221','none'}; d_b = d_b_choice; end
            if jj==4, inversion={'TV_u','norm_S1l1','none'}; d_b = d_b_choice; end
            if jj==5, inversion={'TV_s2','norm_l221','none'}; d_b = d_b_choice; end
            if jj==6, inversion={'TV_s2','norm_S1l1','none'}; d_b = d_b_choice; end
            if jj==7, inversion={'none','norm_l111','CAS8_sym8'}; d_b = d_b_choice; end
            if jj==8, inversion={'TV_u','norm_S1l1','none'}; d_b = []; end

            fprintf('Mask: %s, Testtype: %s, Sim: %d, Parameters setup: %d, Radius: %.2f\n', mask, testtype, sim_string, jj, d_b);

            [I_out, I_acq, mask_out, MR]=wrapper_compressed_acquisition(...
                'im', im_tag, 'ratio', ratio, 'mask', mask,...
                'inv', inversion, 'iter', Nbiter, 'test', testtype,...
                'lambda', lambda_v, 'preproc', preproc, 'tol', tol,...
                'alpha', [], 'sim', sim_string, 'SNR', SNR, 'radius',... 
                d_b, 'idx_metric', 'ssim', 'output', output_folder);

        end
    end
end