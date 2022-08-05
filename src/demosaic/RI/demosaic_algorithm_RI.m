function [rgb_demosaic,time] = demosaic_algorithm_RI(algorithm,I_in,mask,pattern,DynamicRange)

current_folder=pwd;
algorithm_folder=fileparts(mfilename('fullpath'));

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(algorithm, 'RI')
    rgb_mosaic=I_in.*mask;
    cd(fullfile(algorithm_folder,'RI'));
    rgb_demosaic = demosaic_RI(rgb_mosaic,mask,pattern,1,DynamicRange);
    cd(current_folder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(algorithm, 'MLRI')
    rgb_mosaic=I_in.*mask;
    cd(fullfile(algorithm_folder,'MLRI'));
    rgb_demosaic = demosaic_MLRI(rgb_mosaic,mask,pattern,1.4,DynamicRange);
    cd(current_folder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(algorithm, 'MLRI2')
    rgb_mosaic=I_in.*mask;
    cd(fullfile(algorithm_folder,'MLRI2'));
    rgb_demosaic = demosaic_MLRI2(rgb_mosaic,mask,pattern,1,1e-32,DynamicRange);
    cd(current_folder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(algorithm, 'ARI')
    rgb_mosaic=I_in.*mask;
    cd(fullfile(algorithm_folder,'ARI'));
    rgb_demosaic = demosaic_ARI(rgb_mosaic,mask,pattern,DynamicRange);
    cd(current_folder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(algorithm, 'ARI2')
    rgb_mosaic=I_in.*mask;
    cd(fullfile(algorithm_folder,'ARI2'));
    rgb_demosaic = demosaic_ARI2(rgb_mosaic,mask,pattern,DynamicRange);
    cd(current_folder);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=toc;

end
