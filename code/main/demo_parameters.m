function demo_parameters(varargin)

if numel(varargin) >= 1, im_tag = varargin{1}; else, im_tag = 'Beijing_4'; end

% clearvars; 
close all;
fprintf('================ Parametric tests ================\n');

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'jodefu'));

%% Parameters settings - JoDeFu results
% im_tag = 'Beijing_4'; % Image tag
mask = 'BinaryTreeU'; % Mask label
output_folder = 'parameters'; % Output folder
formation_list = 1; % 0 = Spa./spe. degradation, 1 = MRCA, 2 = Mosaicing
test_list = {'lambda', 'tv', 'norm', 'blur'}; % Parameter exploration list

ratio=2; % Scale ratio
Nbiter=250; % Number of algorithm iterations
tol=0; % cost function tolerance (if 0, stop when iteration = Nbiter)

for kk=1:numel(test_list)
    test = test_list{kk};
    output_folder = fullfile(output_folder, test);
    
    preproc='regravg';
    lambda_v = 0.001; % Regression 
    d_b_choice = 1; % Blurring diameter for jodefu v2
    norm_list = {'norm_l221'};
    tv_list = {'TV_c'};
    
    if strcmpi(test, 'lambda'), lambda = logspace(-4,-2.5,15); end
    if strcmpi(test, 'norm'), norm_list = {'norm_l221','norm_S1l1','norm_l211','norm_linf11','norm_l111','norm_Sinfl1'}; end
    if strcmpi(test, 'tv'), tv_list = {'TV_c', 'TV_u', 'TV_s2', 'TV_s3'}; end
    if strcmpi(test, 'blur'), d_b_choice = [1, 1.3:0.1:2]; end

    for jj1 = 1:numel(tv_list)
        for jj2 = 1:numel(norm_list)
            inversion = {tv_list{jj1}, norm_list{jj2}, 'none'};
            for jj3 = 1:numel(d_b_choice)
                d_b = d_b_choice(jj3);
                for ii = formation_list
                    sim_string=0; SNR = [];
                    if ii==0, testtype='nomask'; end
                    if ii==1, testtype='default'; end
                    if ii==2, testtype='msonly'; end
                    if ii==3, testtype='default'; sim_string=1; SNR=25; end

                        fprintf('Mask: %s, Testtype: %s, Sim: %d, Radius: %.2f\n', mask, testtype, sim_string, d_b);

                        [I_out, I_acq, mask_out, MR]=wrapper_compressed_acquisition(...
                            'im', im_tag, 'ratio', ratio, 'mask', mask,...
                            'inv', inversion, 'iter', Nbiter, 'test', testtype,...
                            'lambda', lambda_v, 'preproc', preproc, 'tol', tol,...
                            'alpha', [], 'sim', sim_string, 'SNR', SNR, 'radius',... 
                            d_b, 'idx_metric', 'ssim', 'output', output_folder);
                end
            end
        end
    end
end