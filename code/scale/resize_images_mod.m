%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Generate the low resolution PANchromatic (PAN) and MultiSpectral (MS) images according to the Wald's protocol. 
%           
% Interface:
%           [I_MS_LR, I_PAN_LR] = resize_images_mod(I_MS,I_PAN,ratio,sensor)
%
% Inputs:
%           I_MS:           MS image upsampled at PAN scale;
%           I_PAN:          PAN image;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Integer value;
%           sensor_MS:      String for type of MS sensor (e.g. 'WV2', 'IKONOS').
%           sensor_PAN:     String for type of PAN sensor (e.g. 'WV2', 'IKONOS').
%           Band_selection: Chosen bands.
%
% Outputs:
%           I_MS_LR:        Low Resolution MS image;
%           I_PAN_LR:       Low Resolution PAN image.
% 
% References:
%           [Wald97]    L. Wald, T. Ranchin, and M. Mangolini, “Fusion of satellite images of different spatial resolutions: assessing the quality of resulting images,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 63, no. 6, pp. 691–699, June 1997.
%           [Aiazzi06]  B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591–596, May 2006.
%           [Vivone14]  G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                       IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I_MS_LR, I_PAN_LR] = resize_images_mod(I_MS,I_PAN,ratio,sensor_MS,sensor_PAN,Band_selection)

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'Matlab_toolboxes'));

if nargin<5
    sensor_PAN=sensor_MS;
end

I_MS = double(I_MS);
I_PAN = double(I_PAN);
flag_PAN_MTF = 0;

switch sensor_MS
    case 'QB' 
        flag_resize_new = 2; % MTF usage
        GNyq = [0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
    case 'IKONOS'
        flag_resize_new = 2; % MTF usage
        GNyq = [0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
    case 'GeoEye1' 
        flag_resize_new = 2; % MTF usage
        GNyq = [0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
    case 'WV2'
        flag_resize_new = 2; % MTF usage
        GNyq = [0.35 .* ones(1,7), 0.27];
    case 'ALI'
        flag_resize_new = 2; % MTF usage
        GNyq = [0.29,0.30,0.28,0.29,0.28,0.29,0.25,0.25,0.25];
    case 'HYP'
        % GNyq = [0.27*ones(1,70), 0.29*ones(1,172)];
        %VNIR
        GNyq(1:21)    = 0.27;
        GNyq(22:41)   = 0.28;
        GNyq(42:49)   = 0.26;
        GNyq(50:70)   = 0.26;
        %SWIR
        GNyq(71:100)  = 0.30;
        GNyq(101:130) = 0.30;
        GNyq(131:177) = 0.27;
        GNyq(178:242) = 0.27;
        % Error_Measure = 0.03;
        % GNyq = GNyq + Error_Measure;
        flag_resize_new = 2; % MTF usage
    case 'none'
        if nargin==6
            GNyq = 0.29 .* ones(1,length(Band_selection));
        else
            GNyq = 0.29 .* ones(1,size(I_MS,3));
        end
        flag_resize_new = 1; % Bicubic Interpolator
end
if nargin==6 && ~strcmp(sensor_MS,'none')
    GNyq=GNyq(Band_selection);
end

flag_resize_new_PAN=2;
switch sensor_PAN
    case 'QB'
        GNyqPan = 0.15;
    case 'IKONOS'
        GNyqPan = 0.17;
    case 'GeoEye1'
        GNyqPan = 0.16;
    case 'WV2'
        GNyqPan = 0.11;
    case 'ALI'
        GNyqPan = 0.13;
    case 'none'
        flag_resize_new_PAN =1;
end
if flag_resize_new == 1
    
    %%% Bicubic Interpolator MS
    I_MS_LP = zeros(round(size(I_MS,1)/ratio),round(size(I_MS,2)/ratio),size(I_MS,3));
    
    for idim=1:size(I_MS,3)
        I_MS_LP(:,:,idim) = imresize(I_MS(:,:,idim),1/ratio);
    end
    
    I_MS_LR = double(I_MS_LP);
    
elseif flag_resize_new == 2
    
    %%% MTF
    
    %%% Filtering with sensor MTF MS
    N = 41;
    I_MS_LP = zeros(size(I_MS));
    fcut = 1 / ratio;
    
    for ii = 1 : size(I_MS,3)
        alpha = sqrt(((N-1)*(fcut/2))^2/(-2*log(GNyq(ii))));
        H = fspecial('gaussian', N, alpha);
        Hd = H./max(H(:));
        h = fwind1(Hd,kaiser(N));
        I_MS_LP(:,:,ii) = imfilter(I_MS(:,:,ii),real(h),'replicate');
    end

    %%% Decimation MS
    I_MS_LR = imresize(I_MS_LP,1/ratio,'nearest');
        
end

if flag_resize_new_PAN==2 && flag_PAN_MTF == 1
    %%% Filtering with sensor MTF PAN
    alpha = sqrt(((N-1)*(fcut/2))^2/(-2*log(GNyqPan)));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    I_PAN = imfilter(I_PAN,real(h),'replicate');
    %%% Decimation PAN
    I_PAN_LR = imresize(I_PAN,1/ratio,'nearest');
else
    %%% Bicubic Interpolator PAN
    I_PAN_LR = imresize(I_PAN,1/ratio);
end
    

end

