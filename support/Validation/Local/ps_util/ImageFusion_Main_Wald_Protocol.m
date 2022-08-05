%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  IMAGE FUSION MAIN:  WALD'S PROTOCOL %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%% Analyzed image choice

im_tag = 'WV2'; %sensor = 'none';
sensor = 'WV2';

% im_tag = 'India1'; %sensor = 'none';
% sensor = 'QB';

% im_tag = 'China'; % sensor = 'none';
% sensor = 'IKONOS';

% im_tag = 'Sts1'; 
% sensor = 'none';

% im_tag = 'Sts2'; 
% sensor = 'none';

% im_tag = 'Tls1';
% sensor = 'none';

% im_tag = 'Tls2'; 
% sensor = 'none';


%% Quality Index Blocks
Qblocks_size = 32;

%% Interpolator
bicubic = 0;

%% Cut Final Image
flag_cut_bounds = 1;
dim_cut = 11;

%% Threshold values out of dynamic range
thvalues = 0;

%% Print Eps
print = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%% Dataset load %%%%%%%%%%%%%%%%%%%%%%%%%%
switch im_tag
    case 'WV2'
        load('../data/WVRome.mat');
        I1_loaded = I1;
        I2_loaded = I2;
        im_prepare='resize';
        Resize_fact = 4;
        L = 11; % Radiometric Resolution Sensor
    case 'India1'
        load('Datasets/India1.mat');
        I1_loaded = I1;
        I2_loaded = I2;
        im_prepare='resize';
        Resize_fact = 4;
        L = 11; % Radiometric Resolution Sensor
    case 'China'
        load('../data/IKChina.mat');
        I1 = I1(:,:,[1,3,4,2]);
        I1_loaded = I1;
        I2_loaded = I2;
        im_prepare='resize';
        Resize_fact = 4;
        L = 11; % Radiometric Resolution Sensor
    case 'Sts1'
        I_GT = ReadImage('Datasets/Strsbrg1/Sts1_MS_1024_1m');
        I1 = ReadImage('Datasets/Strsbrg1/Sts1_MS_1024_4m');
        I2 = ReadImage('Datasets/Strsbrg1/Sts1_Pan_1024_1m');
        im_prepare = 'donothing';
        Resize_fact = 4;
        L = 8; % Radiometric Resolution Sensor
    case 'Sts2'
        I_GT = ReadImage('Datasets/Strsbrg2/Sts2_MS_1024_1m');
        I1 = ReadImage('Datasets/Strsbrg2/Sts2_MS_1024_4m');
        I2 = ReadImage('Datasets/Strsbrg2/Sts2_Pan_1024_1m');
        im_prepare = 'donothing';
        Resize_fact = 4;
        L = 8; % Radiometric Resolution Sensor
    case 'Tls1'
        I_GT = ReadImage('Datasets/Tls1/Tls1_MS_1024_1m');
        I1 = ReadImage('Datasets/Tls1/Tls1_MS_1024_4m');
        I2 = ReadImage('Datasets/Tls1/Tls1_Pan_1024_1m');
        im_prepare = 'donothing';
        Resize_fact = 4;
        L = 8; % Radiometric Resolution Sensor
    case 'Tls2'
        I_GT = ReadImage('Datasets/Tls2/Tls2_MS_1024_1m');
        I1 = ReadImage('Datasets/Tls2/Tls2_MS_1024_4m');
        I2 = ReadImage('Datasets/Tls2/Tls2_Pan_1024_1m');
        im_prepare = 'donothing';
        Resize_fact = 4;
        L = 8; % Radiometric Resolution Sensor
end

if strcmp(im_tag,'WV2') || strcmp(im_tag,'India1') || strcmp(im_tag,'China')
    I_GT = double(I1_loaded);
end

%% %%%%%%%%%%%%%    Preparation of image to fuse            %%%%%%%%%%%%%%

if strcmp(im_prepare,'resize')
    if (size(I1_loaded,1) == size(I2_loaded,1)) && (size(I1_loaded,2) == size(I2_loaded,2))
          I1 = resize_images(I1_loaded,I2_loaded,Resize_fact,sensor);
          I2 = double(I2);
    else
          [I1,I2]=resize_images(I1_loaded,I2_loaded,Resize_fact,sensor);
    end
end


%% Upsampling

if strcmp(im_tag,'WV2') || strcmp(im_tag,'India1') || strcmp(im_tag,'China')
    
    I1LR = I1;
    
    if bicubic
        H = zeros(size(I2,1),size(I2,2),size(I1,3));    
        for idim = 1 : size(I1,3)
            H(:,:,idim) = imresize(I1(:,:,idim),Resize_fact);
        end
        I1 = H;
    else
        I1 = interp23tap(I1,Resize_fact);
    end
else
    %%% Sensor Pleiades
    I1LR = imresize(I1,1/Resize_fact,'nearest'); 
end