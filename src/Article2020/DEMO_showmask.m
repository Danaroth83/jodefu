%% Saves Image files

rng('default');  % For reproductible results

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'mosaic'),...
        fullfile(project_folder,'visualization'));
output_folder=fullfile(current_folder,'..','..','data','output','mask_pattern');
    
mask='mindis';

% mask='random';
sizes=[18,18,3];

a=load_mask('mask',mask,'sizes',sizes);
in=visualize_mask(a);

r=9; % magnification size
b=imresize(in,r,'nearest'); 

% Add black borders
b(1:r:end,:,:)=0; 
b(:,1:r:end,:)=0; 
b(r:r:end,:,:)=0; 
b(:,r:r:end,:)=0; 
b=padarray(b,[2,2],0);

mkdir(output_folder);
imshow(b)
imwrite(b,fullfile(output_folder,sprintf('mask_%s.png',mask)));

