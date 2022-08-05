%% MAIN: LOAD DATASETS (REDUCED RESOLUTION)

%
% Author:
% Daniele Picone
%
% Description:
% This function loads a series of multimodal datasets of High and Low 
% resolution (respectively a Panchromatic and a multi/hyperspectral) 
% remotely sensed images representing the same scene
%
% Usage:
% I_req=load_Dataset_Pansharpening_RR(im_tag,'Field',Value);
% Eg:
% - Load dataset with default options
%     I_req=load_Dataset_Pansharpening_RR('Rio_cut2_WV3_WV3_all');
% - Load dataset with simulated PAN from GT data
% I_req=load_Dataset_Pansharpening_RR('Rio_cut2_WV3_WV3_all','flag_PANfromGT',1);
% - Change scale ratio between HR and LR image
%
% in particular it returns a cell I_req containing in the field image:

% MS_LR: The low resolution image (typically multispectral or hyperspectral)
% PAN: The high resolution image (it supposes here a monochromatic acquisition)
% GT: Ground truth (Notice that its scale could be different from I_PAN)
% REF: The reference image for high resolution image

% (Note: If a different subset of images is needed, the list can be changed in the field 'request')

% Each of this images are accompained by a set of metadata whose fields are:

% data: the image itself
% type: the sensor type of the image (panchromatic: 'PAN', multispectral: 'MS', hyperspectral: 'HS')
% request: the requested kind of the image
% ratio: scale ratio with respect to LR image (I_MS_LR)
% IntBits: Bits used to represent the radiometric intensity
% DynamicRange: maximum value the image intensity
% sensor: sensor tag used for the acquisition
% Bands_select: Selection of bands from the original that the user selected to sharpen
% Bands_to_display: Indices of bands for visualization
% GNyq: Gain of the MTF assoiciated to the sensor at the Nyquist frequency;
% KerBlu: The filter used for the image degradation;
% start_pos: The downsampling position of the pixel after degradation
% KerBlu0: A mean of the filters used in KerBlu
% spectralweights: The percentage of overlap of each band at LR on the HR one
% Qblocks_size: Default size for the blocks used for validation;
% edge_cut: Edges' length to be cut on the final image;
% type: Specifies the source of the reference 0=dataset PAN, 1=Multiplatform PAN, 2=Simulated PAN from MS, 3=Multiplatform MS, 4=PAN-ALI, 5= Mean of Multiplatform MS
% flag_resizedoriginal: Specifies if the signal is in its dataset size or was resized
% im_tag: A well formatted im_tag for the image
% GSD: Ground Sample Distance (in meters)
% time: if the image is processed, the time used to process it

% The function takes as input:
% im_tag: A descriptor of the image to load (format:
% [place]_cut[N]_[sensor]_[bands]; ex: SanFrancisco_cut2_ALI_VNIR)

% Other option described below 

function I_req=Load_Dataset_Pansharpening(im_tag,varargin)

disp('Reading dataset...');

current_folder = fileparts(mfilename('fullpath'));
support_folder = fullfile(current_folder,'..');
data_folder    = fullfile(current_folder,'..','..','data','Pansharpening');

addpath(fullfile(support_folder,'Scaling'),...
        fullfile(support_folder,'Filtering'));

%% Note: the variable ratio is the desired ratio between the high and low resolution image after processing, not the actual one

flag_loaddefaultratio=1;
flag_imresize_PAN=1;
flag_imresize_PAN_sim=1;
flag_imresize_MS_RR=0;
flag_imresize_PAN_RR=1;
flag_imresize_HS_RR=0;
flag_imresize_REF=1;
flag_interpolation=0;
flag_avoidPANresize=0;
flag_useREF=0;
flag_histmatchREF=2;
% flag_simulation=0;
flag_cutpower2 = 0;
comp_factor=1;
comp_method='none';
flag_PANfromGT=0;
SNR_PANfromGT=[];
test='RR';
ratio_tgt_PAN=[];
ratio_tgt_REF=[];
request=[];

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if strcmpi(pname,'ratio')                         % Scale ratio between HR and LR image
        ratio=pval;
        if ~isempty(ratio), flag_loaddefaultratio=0; end
    elseif strcmpi(pname,'ratio_tgt_REF')             % Scale ratio between Reference image and LR image
        ratio_tgt_REF=pval;
    elseif strcmpi(pname,'ratio_tgt_PAN')             % Pre-resizing of the PAN signal
        ratio_tgt_PAN=pval;
    elseif strcmpi(pname,'flag_imresize_PAN')         % 1=Prepare PAN image with imresize, 0=with MTF-filtered resizing
        flag_imresize_PAN=pval;
    elseif strcmpi(pname,'flag_imresize_PAN_sim')     % Downsampling of PAN to simulate a different resolution input PAN image
        flag_imresize_PAN_sim=pval;
    elseif strcmpi(pname,'flag_imresize_MS_RR')       % Resizing MS for RR validation 1=with imresize, 0=with MTF-filtered downsampling
        flag_imresize_MS_RR=pval;
    elseif strcmpi(pname,'flag_imresize_PAN_RR')      % Resizing PAN for RR validation 1=with imresize, 0=with MTF-filtered downsampling
        flag_imresize_PAN_RR=pval;
    elseif strcmpi(pname,'flag_imresize_HS_RR')       % Resizing HS for RR validation 1=with imresize, 0=with MTF-filtered downsampling
        flag_imresize_HS_RR=pval;
    elseif strcmpi( pname,'flag_imresize_REF')         % Generation of PAN FR reference image 1=with imresize, 0=with MTF-filtered downsampling
        flag_imresize_REF=pval;
   elseif strcmpi(pname,'flag_interpolation')        % Upsampling option: 1=bicubic (imresize), 0=with a kernel generated with fir1, 2=with a Lagrange polynomial kernel
        flag_interpolation=pval;
    elseif strcmpi(pname,'flag_avoidPANresize')       % If 1, it skips the interpolation of the PAN
        flag_avoidPANresize=pval;
    elseif strcmpi(pname,'flag_useREF')               % Set image as FR reference image 0=dataset PAN, 1=Multiplatform PAN, 2=Simulated PAN from MS, 3=Multiplatform MS, 4=PAN-ALI, 5= Mean of Multiplatform MS 
        flag_useREF=pval;
    elseif strcmpi(pname,'flag_histmatchREF')         % if 1, the reference REF is histogram matched to MS; 2=normalized by radiometric resolution; 3=by max value; 0=non-modified REF 
        flag_histmatchREF=pval;
%     elseif strcmpi(pname,'flag_simulation')         % Use the simulator to create a simulated signal
%         flag_simulation=pval;
    elseif strcmpi(pname,'flag_cutpower2')            % Cut the image sizes to be a power of 2
        flag_cutpower2=pval;
    elseif strcmpi(pname,'comp_factor')               % Radiometrically compresses the image by the factor (must be integer)
        comp_factor=pval;
    elseif strcmpi(pname,'comp_method')               % Chooses how to compress the images before fusion
        comp_method=pval;
    elseif strcmpi(pname,'flag_PANfromGT')            % If 1, simulates the PAN as weighted sum of GT
        flag_PANfromGT=pval;
    elseif strcmpi(pname,'SNR_PANfromGT')             % SNR of the simulated PAN
        SNR_PANfromGT=pval;
    elseif strcmpi(pname,'request')                   % List of requested image to load
        request=pval;
    elseif strcmpi(pname,'test')                      % 'RR'= Reduced Resolution; 'FR': Full Resolution
        test=pval;   
    end
end
if isempty(request)
    if strcmpi(request,'FR'), request={'MS_LR','PAN','REF'}; 
    else, request={'MS_LR','PAN','GT'}; end
end
if ~iscell(request), request={request}; end

% im_prepare='resize';   % when imprepare='donothing' it doesn't generate MS-LR by degrading GT, but it is loaded from .mat file

if flag_useREF==2 || flag_useREF==3 || flag_useREF==5
    flag_useREF=1;
    if exist('flag_RR_val','var') && flag_RR_val==5
        disp('Warning: This script doesn''t allow MS as reference');
    end
end

sensor_HS_list={'HYP','AVIRIS','CHR','ROSIS'};

[place,cut,Bands,sensor,sensor_PAN]=parse_im_tag(im_tag);
% type='HS';
[vcut,hcut,edge_cut,edge_cut_PAN,Qblocks_size]=load_cut(place,cut,sensor,sensor_PAN,'PAN');
[Bands_to_sharpen,Bands_to_display]=load_Bands_to_sharpen(sensor,Bands,place);
[~,Bands_to_display]=min(abs(repmat(Bands_to_sharpen',[1,3])-repmat(Bands_to_display,[length(Bands_to_sharpen),1])),[],1);
[GSD_GT,L]=load_resolution( sensor, im_tag, 'MS');
ratio_PAN=GSD_GT/load_resolution(sensor_PAN,im_tag,'PAN');
Band_overlap_PAN=load_Band_overlap(sensor,sensor_PAN,'HS','PAN');
Band_overlap_PAN=unique(find(ismember(Bands_to_sharpen,Band_overlap_PAN)));
if flag_loaddefaultratio==1, ratio=ratio_PAN; end


%% Dataset options
if strncmpi(place,'Rome',4) || strncmpi(place,'China',5) || strncmpi(place,'Pleiades',8)
    % Those are just available at RR as the PAN is resized
    if flag_loaddefaultratio==1, ratio=4; end   
    filename=[place,'_RR'];
    load(fullfile(data_folder,[filename,'.mat']),'I_GT','I_PAN');
    [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN,vcut,hcut,1,edge_cut_PAN,edge_cut_PAN);
    [I_MS_loaded,edge_cut]=imcrop_custom(I_GT(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
    
elseif strncmpi(place,'CHRIS',5) || strncmpi(place,'Paris',5)
    

    load(fullfile(data_folder,[place,'.mat']),'I_MS_LR','I_PAN');
    [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
    [I_MS_loaded,edge_cut]=imcrop_custom(I_MS_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);

elseif strncmpi(place,'Pavia',5) || strncmpi(place,'Moffett',7)
    % Hyperspectral Datasets. PAN is not available
    
    if flag_loaddefaultratio==1, ratio=5; end
    filename=fullfile('..','HS',place);
    load(fullfile(data_folder,[filename,'.mat']),'I_MS');
    [I_MS_loaded,edge_cut]=imcrop_custom(I_MS(:,:,Bands_to_sharpen),vcut,hcut,1,0,0);
    flag_PANfromGT=1;

elseif strncmpi(place,'Washington',10) || strncmpi(place,'Janeiro',7) || ...
        strncmpi(place,'Stockholm',9) || strncmpi(place,'Tripoli',7) || ...
        strncmpi(place,'Vancouver',9) || strncmpi(place,'Hobart',6)
    % Various datasets from DigitalGlobe
    load(fullfile(data_folder,[place,'.mat']),'I_MS_LR','I_PAN');
    [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
    [I_MS_loaded,edge_cut]=imcrop_custom(I_MS_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
    
elseif strncmpi(place,'Toulouse',8) || strncmpi(place,'Indianapolis',12)
    % Datasets with double cut available
    filename=fullfile(data_folder,[place,'.mat']);
    if strcmpi(cut,'cut2') % Hack to load a second image
        load(filename,'I_PAN2','I_MS2_LR'); I_PAN=I_PAN2; I_MS_LR=I_MS2_LR; clear I_PAN2; clear I_MS2_LR;
    else
        load(filename,'I_PAN','I_MS_LR');
    end
    [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
    [I_MS_loaded,edge_cut]=imcrop_custom(I_MS_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);

elseif strncmpi(place,'Sofia',5) || strncmpi(place,'Sudbury',7) || strncmpi(place,'Fuji',4) 
    % Hyperion ALI datasets. Other than the PAN, we have both the MS (ALI) and the HS (Hyperion)   
    load_file=fullfile(data_folder,[place,'.mat']);
    if strcmpi(sensor,'ALI') && strcmpi(sensor_PAN,'ALI')
        load(load_file,'I_MS0_LR','I_PAN0');
        [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN0,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
        [I_MS_loaded,edge_cut]=imcrop_custom(I_MS0_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
    elseif strcmpi(sensor,'HYP') && strcmpi(sensor_PAN,'ALI')
        load(load_file,'I_MS_LR','I_PAN');
        [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
        [I_MS_loaded,edge_cut]=imcrop_custom(I_MS_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
    end
    
elseif strncmpi(place,'Beijing',7) || strncmpi(place,'RdJ',3) || ...
        strncmpi(place,'Sydney',6) || strncmpi(place,'SaoPaulo',8) || ...
        strncmpi(place,'SanFrancisco',12)
    % Multiplatform Datasets, taken from Hyperion ALI and various commercial sensors
    
    if flag_loaddefaultratio==1 && ratio>4, ratio=3; end %% Hack for multiplatform ratios
    
    load_file=fullfile(data_folder,'..','Multiplatform',[place,'.mat']);
    if ~strcmpi(sensor,'HYP') && ~strcmpi(sensor_PAN,'ALI') && ~strcmpi(sensor,'ALI')
        load(load_file,'I_MS0_LR','I_PAN0');
        [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN0,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
        [I_MS_loaded,edge_cut]=imcrop_custom(I_MS0_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
        flag_useREF=0;
    elseif strcmpi(sensor,'HYP') && ~strcmpi(sensor_PAN,'ALI')
        load(load_file,'I_MS2_LR','I_PAN');
        [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
        [I_MS_loaded,edge_cut]=imcrop_custom(I_MS2_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
    elseif strcmpi(sensor,'ALI') && strcmpi(sensor_PAN,'ALI')
        load(load_file,'I_MS1_LR','I_PAN1');
        [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN1,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
        [I_MS_loaded,edge_cut]=imcrop_custom(I_MS1_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
        flag_useREF=0;
    elseif strcmpi(sensor,'HYP') && strcmpi(sensor_PAN,'ALI')
        load(load_file,'I_MS2_LR','I_PAN2');
        [I_PAN_loaded,edge_cut_PAN]=imcrop_custom(I_PAN2,vcut,hcut,ratio_PAN,edge_cut_PAN,edge_cut_PAN);
        [I_MS_loaded,edge_cut]=imcrop_custom(I_MS2_LR(:,:,Bands_to_sharpen),vcut,hcut,1,edge_cut,edge_cut);
    end
    
    if flag_useREF==4
        sensor_REF='ALI';
        load(load_file,'I_PAN2');
        I_REF=I_PAN2;
    else
        [~,sensor_REF]=load_default_sensor(place);
        load(load_file,'I_PAN');
        I_REF=I_PAN;
    end
    [~,~,~,edge_cut_REF,~]=load_cut(im_tag,cut,sensor,sensor_REF,'PAN');
    Band_overlap_REF=load_Band_overlap(sensor,sensor_REF,'MS','PAN');
    [~,L_REF]=load_resolution(sensor_REF,place);
    ratio_REF=load_resolution(sensor,im_tag)/load_resolution(sensor_REF,im_tag,'PAN');
    [I_REF,edge_cut_REF]=imcrop_custom(I_REF,vcut,hcut,ratio_REF,edge_cut_REF,edge_cut_REF);
    I_REF = double(I_REF);
    
end

%% Sensor fix (To remove in later version; it sets sensors to 'none' for non-implemented MTF filters)
if isempty(Bands_to_sharpen) || any(isnan(Bands_to_sharpen)), Bands_to_sharpen=1:size(I_MS_loaded,3); end
if isempty(Band_overlap_PAN) || any(isnan(Band_overlap_PAN)), Band_overlap_PAN=1:size(I_MS_loaded,3); end
if isempty(Bands_to_display) || any(isnan(Bands_to_display)), Bands_to_display=1:3; end

%% Loading MTF gains at Nyquist frequencies
GNyq_PAN=load_MTF_PAN(sensor_PAN);
GNyq=load_MTF(sensor, im_tag, Bands_to_sharpen);

%% PAN resizing
if isempty(ratio_tgt_REF), ratio_tgt_REF=ratio; end
if isempty(ratio_tgt_PAN), ratio_tgt_PAN=ratio; end
GNyq2_PAN=mean(load_MTF(sensor,im_tag,intersect(1:length(Bands_to_sharpen),Band_overlap_PAN)));
if flag_PANfromGT==1
    if isnan(load_resolution(sensor,im_tag,'PAN')), sensor_PAN='none'; else, sensor_PAN=sensor; end
    ratio_PAN=ratio_tgt_PAN/ratio_tgt_REF*ratio;
    spectralweights=load_spectralweights(sensor,sensor_PAN,Bands_to_sharpen,'MS');
    ratio_tgt_PAN=min(ratio_tgt_PAN,ratio_tgt_REF);
    [I_temp,edge_cut_PAN]=scale_image(I_MS_loaded,ratio_tgt_PAN/ratio_tgt_REF,flag_interpolation,...
        flag_imresize_REF,'GNyq',GNyq,'edge',edge_cut,'removeedge',0);
    I_PAN_loaded=permute(sum(double(I_temp).*shiftdim(spectralweights,-2),3),[1,2,4,3]);
    clear I_temp;
    if ~isempty(SNR_PANfromGT)
        SNR_nn=10^(SNR_PANfromGT/10);
        Var_Noise=sum(I_PAN_loaded(:).^2)/numel(I_PAN_loaded)/SNR_nn;
        Random_Noise=sqrt(Var_Noise)*randn(size(I_PAN_loaded));
        I_PAN_loaded=I_PAN_loaded+Random_Noise;
        I_PAN_loaded(I_PAN_loaded<=0)=0;
        I_PAN_loaded(I_PAN_loaded>=2^L-1)=2^L-1;
        % I_PAN_loaded=I_PAN_loaded+randn(size(I_PAN_loaded))*sqrt(mean((I_PAN_loaded(:)-mean(I_PAN_loaded(:))).^2)/(10^(SNR_PANfromGT/10)));
    end
end

if flag_avoidPANresize~=1
    [I_PAN_loaded,edge_cut_PAN]=scale_image(I_PAN_loaded,ratio/ratio_tgt_PAN,flag_interpolation,...
        flag_imresize_PAN,'GNyq',GNyq_PAN,'GNyq2',GNyq2_PAN,...
        'edge',edge_cut_PAN,'removeedge',0);
    if ratio_tgt_PAN~=ratio, flag_resizedoriginal_PAN=1; else, flag_resizedoriginal_PAN=0; end
else
    I_PAN_loaded=I_PAN_loaded(edge_cut_PAN+1:end-edge_cut_PAN,edge_cut_PAN+1:end-edge_cut_PAN);
    flag_resizedoriginal_PAN=0;
end

%% Change variable type
I_MS_loaded=double(I_MS_loaded);
I_PAN_loaded=double(I_PAN_loaded);

%% Fix for missing reference
if ~exist('sensor_REF','var') || flag_useREF==0
    flag_useREF=0;
    sensor_REF=sensor_PAN;
    ratio_REF=ratio_PAN;
    edge_cut_REF=edge_cut_PAN;
    Band_overlap_REF=Band_overlap_PAN;
    L_REF=L;
    I_REF=I_PAN_loaded;
end

%% Ground Truth, MS and PAN image RR construction
ratio_actual_PAN=(size(I_PAN_loaded,1)-2*edge_cut_PAN)/(size(I_MS_loaded,1)-2*edge_cut);
ratio_actual_REF=(size(I_REF,1)-2*edge_cut_REF)/(size(I_MS_loaded,1)-2*edge_cut);
if any(strcmpi(sensor,sensor_HS_list)), flag_imresize_RR=flag_imresize_HS_RR; else, flag_imresize_RR=flag_imresize_MS_RR; end

if strcmpi(test,'FR')
    [I_MS_loaded,I_PAN_loaded,I_REF,edge_temp]=size_check(...
        [1,ratio_actual_PAN/ratio,ratio_actual_REF/ratio_tgt_REF],...
        I_MS_loaded,I_PAN_loaded,I_REF,'edge',[edge_cut,edge_cut_PAN,edge_cut_REF]);
    edge_cut=edge_temp(1); edge_cut_PAN=edge_temp(2); edge_cut_REF=edge_temp(3);
    I_GT=I_MS_loaded(edge_cut+1:end-edge_cut,edge_cut+1:end-edge_cut,:);
    I_MS_LR_e=I_MS_loaded;
else
    [I_MS_loaded,I_PAN_loaded,I_REF,edge_temp]=size_check(...
        [ratio_tgt_REF,ratio_tgt_REF*ratio_actual_PAN/ratio,ratio_actual_REF],...
        I_MS_loaded,I_PAN_loaded,I_REF,'edge',[edge_cut,edge_cut_PAN,edge_cut_REF]);
    edge_cut=edge_temp(1); edge_cut_PAN=edge_temp(2); edge_cut_REF=edge_temp(3);
    I_GT=I_MS_loaded(edge_cut+1:end-edge_cut,edge_cut+1:end-edge_cut,:);
    [I_MS_LR_e,edge_cut]=scale_image(I_MS_loaded,1/ratio_tgt_REF,flag_interpolation,...
        flag_imresize_RR,'GNyq',load_MTF(sensor,im_tag,Bands_to_sharpen),...
        'edge',edge_cut,'removeedge',0);
end

% PAN resizing

if strcmpi(test,'FR')
    [I_PAN,edge_cut_PAN]=scale_image(I_PAN_loaded,ratio_tgt_PAN/ratio_actual_PAN,...
        flag_interpolation,flag_imresize_PAN_sim,'GNyq',GNyq_PAN,'GNyq2',...
        GNyq2_PAN,'edge',edge_cut_PAN,'removeedge',0);
    I_PAN=scale_image(I_PAN,ratio/ratio_tgt_PAN,...
        flag_interpolation,flag_imresize_PAN_sim,'GNyq',GNyq_PAN,'GNyq2',...
        GNyq2_PAN,'edge',edge_cut_PAN,'removeedge',1);
else
    ratio1_tgt_PAN=min(ratio_tgt_PAN,ratio_actual_PAN);
    [I_PAN,edge_cut_PAN]=scale_image(I_PAN_loaded,ratio1_tgt_PAN/ratio_actual_PAN,...
        flag_interpolation,flag_imresize_PAN_sim,'GNyq',GNyq_PAN,'GNyq2',...
        GNyq2_PAN,'edge',edge_cut_PAN,'removeedge',0);
    [I_PAN,edge_cut_PAN]=scale_image(I_PAN,...
    (ratio_tgt_PAN/ratio1_tgt_PAN)/ratio_tgt_REF,flag_interpolation,...
    flag_imresize_PAN_RR,'GNyq',GNyq_PAN,'GNyq2',GNyq2_PAN,...
    'edge',edge_cut_PAN,'removeedge',0);
    I_PAN=scale_image(I_PAN,ratio/ratio_tgt_PAN,flag_interpolation,...
    flag_imresize_PAN_RR,'GNyq',GNyq_PAN,'GNyq2',GNyq2_PAN,...
    'edge',edge_cut_PAN,'removeedge',1);
end
edge_cut_PAN=0;

%% FR Reference resize
if strcmpi(test,'FR')
    I_REF=scale_image(I_REF,ratio_tgt_REF/ratio_actual_REF,flag_interpolation,flag_imresize_REF,...
        'GNyq',load_MTF_PAN(sensor_REF),'edge',edge_cut_REF,'removeedge',1);
else
    I_REF=scale_image(I_REF,1/ratio_actual_REF,flag_interpolation,flag_imresize_REF,...
        'GNyq',load_MTF_PAN(sensor_REF),'edge',edge_cut_REF,'removeedge',1);
end

%% Compression
if flag_PANfromGT==1,L_PAN=L; else, [~,L_PAN]=load_resolution(sensor_PAN,im_tag,'PAN'); end
I_PAN=compress_image(I_PAN,comp_factor,comp_method,L_PAN);
I_MS_LR_e=compress_image(I_MS_LR_e,comp_factor,comp_method,L);

%% Cut images
if flag_cutpower2==1
    [L1,L2,~]=size(I_MS_LR);
    dim1=2^(nextpow2(L1-1)-1);
    dim2=2^(nextpow2(L2-1)-1);
    vcut=1:dim1; % hcut=(1:dim2)+floor((size(I_MS_LR,1)-dim1)/2);
    hcut=1:dim2; % vcut=(1:dim1)+floor((size(I_MS_LR,2)-dim2)/2);
    [I_MS_LR_e,edge_cut,~]=imcrop_custom(I_MS_LR_e,vcut,hcut,1,edge_cut,edge_cut);
    [~,~,I_PAN]=imcrop_custom(I_PAN,vcut,hcut,ratio,0,0);
    [~,~,I_REF]=imcrop_custom(I_REF,vcut,hcut,ratio_tgt_REF,0,0);
    [~,~,I_GT]=imcrop_custom(I_GT,vcut,hcut,ratio_tgt_REF,0,0);
    if L1~=dim1, fprintf('Vertical GT size changed to %d, so LR image size is a power of 2\n',dim1*ratio_tgt_REF); end
    if L2~=dim2, fprintf('Horizontal GT size changed to %d, so LR image size is a power of 2\n',dim2*ratio_tgt_REF); end
end

%% FR Reference histogram matching

if flag_histmatchREF==1
    I_EXP=scale_image(I_MS_LR_e,ratio,flag_interpolation,[],'edge',edge_cut);
    onemat=ones(size(I_REF));
    I_MS_temp=I_EXP(:,:,ismember(Bands_to_sharpen, intersect(Bands_to_sharpen,Band_overlap_REF)));
    I_REF=(I_REF-mean2(I_REF)*onemat)/std2(I_REF)*std2(I_MS_temp)+mean2(I_MS_temp)*onemat;
    clear I_MS_temp
elseif flag_histmatchREF==2
    I_REF=I_REF*2^L/2^L_REF;
elseif flag_histmatchREF==3
    I_REF=I_REF*max(I_MS(:))/max(I_REF(:));
end

%% Final metadata
if length(Bands_to_sharpen)~=2^(nextpow2(length(Bands_to_sharpen)))
    disp('Q2n index is advised to work with 2^n bands (where n is an integer)');
end
if any(strcmp(sensor_HS_list,sensor)), type='HS'; else, type='MS'; end
[KerBlu,start_pos]=load_KerBlu(im_tag,ratio,GNyq,min(size(I_PAN,1),size(I_PAN,2)),false);
[KerBlu0,~]=load_KerBlu(im_tag,ratio,GNyq,min(size(I_PAN,1),size(I_PAN,2)));
spectralweights=load_spectralweights(sensor,sensor_PAN,Bands_to_sharpen,'MS');
[wavelength_central,wavelength_bw]=load_wavelength( sensor,type,Bands_to_sharpen );
[wavelength_central_PAN,wavelength_bw_PAN]=load_wavelength( sensor_PAN,'PAN' );

%% Assignment of loaded images to the output

n_req=numel(request);
I_req=cell(1,n_req);

for ii=1:n_req

    switch request{ii}

        case {'MS_LR','MS-LR','LR'}
            
            I_MS_LR=I_MS_LR_e(edge_cut+1:end-edge_cut,edge_cut+1:end-edge_cut,:);
            
            I_req{ii}.data=I_MS_LR;
            I_req{ii}.type=type;
            I_req{ii}.label='MS_LR';
            I_req{ii}.IntBits=L;
            I_req{ii}.DynamicRange=2^L-1;
            I_req{ii}.sensor=sensor;
            I_req{ii}.ratio=1;
            I_req{ii}.Bands_to_display=Bands_to_display;
            I_req{ii}.Bands_selected=Bands_to_sharpen;
            I_req{ii}.Band_overlap_PAN=Band_overlap_PAN;
            I_req{ii}.GNyq=GNyq;
            I_req{ii}.Qblocks_size=Qblocks_size;
            I_req{ii}.edge_cut=0;
            I_req{ii}.im_tag=im_tag;
            I_req{ii}.comp_factor=comp_factor;
            I_req{ii}.comp_method=comp_method;
            I_req{ii}.KerBlu=KerBlu;
            I_req{ii}.KerBlu0=KerBlu0;
            I_req{ii}.start_pos=start_pos;
            I_req{ii}.spectralweights=spectralweights;
            if strcmpi(test,'FR'), I_req{ii}.GSD=GSD_GT; else, I_req{ii}.GSD=GSD_GT*ratio_tgt_REF; end
            I_req{ii}.wavelength=wavelength_central;
            I_req{ii}.bandwidth=wavelength_bw;
            I_req{ii}.size=size(I_req{ii}.data); I_req{ii}.size(3)=size(I_req{ii}.data,3);
            
        case {'MS_LR_e','MS-LR_e'}

            I_req{ii}.data=I_MS_LR_e;
            I_req{ii}.type=type;
            I_req{ii}.label='MS_LR_e';
            I_req{ii}.IntBits=L;
            I_req{ii}.DynamicRange=2^L-1;
            I_req{ii}.sensor=sensor;
            I_req{ii}.ratio=1;
            I_req{ii}.Bands_to_display=Bands_to_display;
            I_req{ii}.Bands_selected=Bands_to_sharpen;
            I_req{ii}.Band_overlap_PAN=Band_overlap_PAN;
            I_req{ii}.GNyq=GNyq;
            I_req{ii}.Qblocks_size=Qblocks_size;
            I_req{ii}.edge_cut=edge_cut;
            I_req{ii}.im_tag=im_tag;
            I_req{ii}.comp_factor=comp_factor;
            I_req{ii}.comp_method=comp_method;
            I_req{ii}.KerBlu=KerBlu;
            I_req{ii}.KerBlu0=KerBlu0;
            I_req{ii}.start_pos=start_pos;
            I_req{ii}.spectralweights=spectralweights;
            if strcmpi(test,'FR'), I_req{ii}.GSD=GSD_GT*ratio_tgt_REF; else, I_req{ii}.GSD=GSD_GT; end
            I_req{ii}.wavelength=wavelength_central;
            I_req{ii}.bandwidth=wavelength_bw;
            I_req{ii}.size=size(I_req{ii}.data); I_req{ii}.size(3)=size(I_req{ii}.data,3);
            
        case 'PAN'
            
            I_req{ii}.data=I_PAN;
            I_req{ii}.type='PAN';
            I_req{ii}.label='PAN';
            I_req{ii}.IntBits=L_PAN;
            I_req{ii}.DynamicRange=2^L_PAN-1;
            I_req{ii}.sensor=sensor_PAN;
            I_req{ii}.Bands_to_display=[1,1,1];
            I_req{ii}.ratio=ratio;
            I_req{ii}.flag_resizedoriginal=flag_resizedoriginal_PAN;
            I_req{ii}.flag_PANfromGT=flag_PANfromGT;
            I_req{ii}.SNR_PANfromGT=SNR_PANfromGT;
            I_req{ii}.edge_cut=edge_cut_PAN;
            I_req{ii}.im_tag=im_tag;
            I_req{ii}.comp_factor=comp_factor;
            I_req{ii}.comp_method=comp_method;
            I_req{ii}.spectralweights=spectralweights;
            if strcmpi(test,'FR'), I_req{ii}.GSD=GSD_GT/ratio; else, I_req{ii}.GSD=GSD_GT*ratio_tgt_REF/ratio; end
            I_req{ii}.wavelength=wavelength_central_PAN;
            I_req{ii}.bandwidth=wavelength_bw_PAN;
            I_req{ii}.size=size(I_req{ii}.data); I_req{ii}.size(3)=size(I_req{ii}.data,3);

        case 'REF'
            
            [wavelength_central_REF,wavelength_bw_REF]=load_wavelength( sensor_REF,'PAN' );
            
            I_req{ii}.data=I_REF;
            I_req{ii}.type='PAN';
            I_req{ii}.label='REF';
            I_req{ii}.IntBits=L_REF;
            I_req{ii}.DynamicRange=2^L_REF-1;
            I_req{ii}.sensor=sensor_REF;
            I_req{ii}.Bands_to_display=[1,1,1];
            I_req{ii}.ratio=ratio_REF;
            I_req{ii}.type=flag_useREF;
            I_req{ii}.flag_histmatchREF=flag_histmatchREF;
            I_req{ii}.im_tag=im_tag;
            I_req{ii}.edge_cut=0;
            if strcmpi(test,'FR'), I_req{ii}.GSD=GSD_GT/ratio_tgt_REF; else, I_req{ii}.GSD=GSD_GT; end
            I_req{ii}.wavelength=wavelength_central_REF;
            I_req{ii}.bandwidth=wavelength_bw_REF;
            I_req{ii}.size=size(I_req{ii}.data); I_req{ii}.size(3)=size(I_req{ii}.data,3);
            
        case 'GT'
            
            I_req{ii}.data=I_GT;
            I_req{ii}.type=type;
            I_req{ii}.label='GT';
            I_req{ii}.IntBits=L;
            I_req{ii}.DynamicRange=2^L-1;
            I_req{ii}.sensor=sensor;
            if strcmpi(test,'FR'), I_req{ii}.ratio=1; else, I_req{ii}.ratio=ratio_tgt_REF; end
            I_req{ii}.Bands_to_display=Bands_to_display;
            I_req{ii}.Bands_selected=Bands_to_sharpen;
            I_req{ii}.Band_overlap_PAN=Band_overlap_PAN;
            I_req{ii}.GNyq=GNyq;
            I_req{ii}.Qblocks_size=Qblocks_size;
            I_req{ii}.im_tag=im_tag;
            I_req{ii}.edge_cut=0;
            I_req{ii}.KerBlu=KerBlu;
            I_req{ii}.KerBlu0=KerBlu0;
            I_req{ii}.start_pos=start_pos;
            I_req{ii}.spectralweights=spectralweights;
            I_req{ii}.GSD=GSD_GT;
            I_req{ii}.wavelength=wavelength_central;
            I_req{ii}.bandwidth=wavelength_bw;
            I_req{ii}.size=size(I_req{ii}.data); I_req{ii}.size(3)=size(I_req{ii}.data,3);

        case {'EXP','MS'}
            
            [I_EXP,~,time_EXP]=scale_image(I_MS_LR_e,ratio,flag_interpolation,...
                [],'edge',edge_cut);
            
            I_req{ii}.data=I_EXP;
            I_req{ii}.type=type;
            I_req{ii}.label='EXP';
            I_req{ii}.IntBits=L;
            I_req{ii}.DynamicRange=2^L-1;
            I_req{ii}.sensor=sensor;
            I_req{ii}.Bands_to_display=Bands_to_display;
            I_req{ii}.Bands_selected=Bands_to_sharpen;
            I_req{ii}.Band_overlap_PAN=Band_overlap_PAN;
            I_req{ii}.ratio=ratio;
            I_req{ii}.im_tag=im_tag;
            I_req{ii}.edge_cut=0;
            I_req{ii}.time=time_EXP;
            I_req{ii}.GNyq=GNyq;
            I_req{ii}.KerBlu=KerBlu;
            I_req{ii}.KerBlu0=KerBlu0;
            I_req{ii}.start_pos=start_pos;
            I_req{ii}.spectralweights=spectralweights;
            if strcmpi(test,'FR'), I_req{ii}.GSD=GSD_GT/ratio; else, I_req{ii}.GSD=GSD_GT*ratio_tgt_REF/ratio; end
            I_req{ii}.wavelength=wavelength_central;
            I_req{ii}.bandwidth=wavelength_bw;
            I_req{ii}.size=size(I_req{ii}.data); I_req{ii}.size(3)=size(I_req{ii}.data,3);
            I_req{ii}.interpolation=flag_interpolation;

    end
        
end

disp('Done');
