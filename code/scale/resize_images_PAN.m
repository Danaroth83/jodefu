%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Generate the low resolution PANchromatic (PAN) images according to the Wald's protocol. 
%           
% Interface:
%           [I_PAN_LR] = resize_images_MS(I_PAN,ratio,GNyq,flag_imresize)
%
% Inputs:
%           I_PAN:          PAN image;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Integer value;
%           GNyqPan:        MTF gain at Nyquist frequency;
%           flag_MTFresize: Resizing method (0,=imresize, 1=MTF-aided)
%
% Outputs:
%           I_PAN_LR:       Low Resolution PAN image;
% 
% References:
%           [Wald97]    L. Wald, T. Ranchin, and M. Mangolini, “Fusion of satellite images of different spatial resolutions: assessing the quality of resulting images,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 63, no. 6, pp. 691–699, June 1997.
%           [Aiazzi06]  B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591–596, May 2006.
%           [Vivone14]  G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                       IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I_PAN_LR] = resize_images_PAN(I_PAN,ratio,GNyqPan,flag_imresize, flag_decimate)

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'Matlab_toolboxes'));

if nargin<=4, flag_decimate=1; end
if nargin<=3, flag_imresize=1; end
if nargin<=2 || isempty(GNyqPan), GNyqPan=0.15; end

I_PAN = double(I_PAN);

if isequal(flag_imresize,1) || strcmpi(flag_imresize,'imresize')
    %%% Bicubic Interpolator PAN
    I_PAN_LR = imresize(I_PAN,1/ratio);
    if flag_decimate==0
        I_PAN_LR=imresize(I_PAN_LR,ratio);
    end
else
    N = 41;
    fcut = 1 / ratio;
    %%% Filtering with sensor MTF PAN
    alpha = sqrt(((N-1)*(fcut/2))^2/(-2*log(GNyqPan)));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    I_PAN_LR = imfilter(I_PAN,real(h),'replicate');
    %%% Decimation PAN
    if flag_decimate==1
        I_PAN_LR = imresize(I_PAN_LR,1/ratio,'nearest');
    end
end