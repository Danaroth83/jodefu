close all; clearvars;

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'jodefu'));

%% Image formation - Software compression

im_tag = 'Washington_4'; % Image tag
compression_ratio = 0.25; % Compression ratio
output_folder = 'formation_software'; % Output folder
methods_list = {'BIN', 'JPEG'}; % Software compression methods

wrapper_software_compression('im_tag', im_tag, 'compression_ratio', ...
    compression_ratio, 'methods', methods_list, 'flag_vis', 1, ...
    'output', output_folder);