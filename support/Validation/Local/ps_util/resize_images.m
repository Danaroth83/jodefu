%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Generate the low resolution PANchromatic (PAN) and MultiSpectral (MS) images according to the Wald's protocol. 
%           
% Interface:
%           [I_MS_LR, I_PAN_LR] = resize_images(I_MS,I_PAN,ratio,sensor)
%
% Inputs:
%           I_MS:           MS image upsampled at PAN scale;
%           I_PAN:          PAN image;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Integer value;
%           sensor:         String for type of sensor (e.g. 'WV2', 'IKONOS').
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
function [I_MS_LR, I_PAN_LR] = resize_images(I_MS,I_PAN,ratio,sensor,varargin)

I_MS = double(I_MS);
I_PAN = double(I_PAN);
flag_PAN_MTF = 0;
flag_resize_new=2;

if strcmp(sensor,'none')
        flag_resize_new = 1; % Bicubic Interpolator
end

if flag_resize_new == 1
    
    %%% Bicubic Interpolator MS
    I_MS_LP = zeros(round(size(I_MS,1)/ratio),round(size(I_MS,2)/ratio),size(I_MS,3));
    
    for idim=1:size(I_MS,3)
        I_MS_LP(:,:,idim) = imresize(I_MS(:,:,idim),1/ratio);
    end
    
    I_MS_LR = double(I_MS_LP);
    
    %%% Bicubic Interpolator PAN
    I_PAN_LR = imresize(I_PAN,1/ratio);
    
elseif flag_resize_new == 2
    
    %%% MTF
    I_MS_LP = MTF(I_MS,sensor,[],ratio);
    %%% Decimation MS
    I_MS_LR = imresize(I_MS_LP,1/ratio,'nearest');
    
    if flag_PAN_MTF == 1
        I_PAN = MTF_PAN(I_PAN,sensor,ratio);
        %%% Decimation PAN
        I_PAN_LR = imresize(I_PAN,1/ratio,'nearest');
    else
        I_PAN_LR = imresize(I_PAN,1/ratio);
    end

        
end

end