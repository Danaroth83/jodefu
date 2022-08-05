close all; clearvars;

current_folder=fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder,'..','jodefu'));

im_tag = 'Washington_4';
compression_ratio = 0.25;
output_folder = 'formation_software';
methods_list = {'BIN', 'JPEG'};

wrapper_software_compression('im_tag', im_tag, 'compression_ratio', ...
    compression_ratio, 'methods', methods_list, 'flag_vis', 1, ...
    'output', output_folder);