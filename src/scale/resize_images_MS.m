%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Generate the low resolution MultiSpectral (MS) images according to the Wald's protocol. 
%           
% Interface:
%           [I_MS_LR] = resize_images_MS(I_MS,ratio,GNyq,flag_resize_new)
%
% Inputs:
%           I_MS:            MS image upsampled at PAN scale;
%           ratio:           Scale ratio between MS and PAN. Pre-condition: Integer value;
%                            It can be in the format [vertical,horizontal]
%           GNyq:            MTF gain at Nyquist frequency;
%           flag_MTFresize:  Resizing method (0,=imresize, 1=MTF-aided)
%
% Outputs:
%           I_MS_LR:        Low Resolution MS image;
% 
% References:
%           [Wald97]    L. Wald, T. Ranchin, and M. Mangolini, “Fusion of satellite images of different spatial resolutions: assessing the quality of resulting images,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 63, no. 6, pp. 691–699, June 1997.
%           [Aiazzi06]  B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591–596, May 2006.
%           [Vivone14]  G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                       IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_MS_LR] = resize_images_MS(I_MS,ratio,GNyq,flag_imresize,flag_decimate,flag_disalign)

current_folder=fileparts(mfilename('fullpath'));
project_folder=fullfile(current_folder,'..');
addpath(fullfile(project_folder,'Matlab_toolboxes'));

if numel(ratio)==1, ratio=[ratio,ratio]; end
if nargin<=2 || isempty(GNyq), GNyq=0.29 .* ones(1,size(I_MS,3)); end
if nargin<=3, flag_imresize=0; end
if nargin<=4, flag_decimate=1; end
if nargin<=5, flag_disalign=1; end

I_MS = double(I_MS);
[L1,L2,Nb]=size(I_MS);

% Fix for even ratios (only works with integer ratios)
if flag_decimate==1 && flag_disalign==0
    interp_fix=2-rem(ratio,2);
    if any(rem(interp_fix,1)~=0)
        disp('Warning: downsampling with non integer scale ratios causes disalignments');
    end
    interp_fix(interp_fix~=2)=1;
    if ~any(rem(interp_fix,1)~=0)
        I_MS=interp23tap_offcenter(I_MS,interp_fix);
        I_MS=I_MS(1:interp_fix(1):end,1:interp_fix(2):end,:);
    end
end

if isequal(flag_imresize,1) || strcmpi(flag_imresize,'imresize')
    
    %%% Bicubic Interpolator MS
    
    I_MS_LR = imresize(I_MS,round([L1,L2]./ratio));
    
    I_MS_LR = double(I_MS_LR);
    if flag_decimate==0
        I_MS_LR=imresize(I_MS_LR,[L1,L2]);
    end
    
else
    
    %%% MTF
    
    %%% Filtering with sensor MTF MS
    N = 41;
    I_MS_LR = zeros([L1,L2,Nb]);
    fcut = 1 ./ ratio;
    
    for ii = 1 : Nb
        alpha = sqrt(((N-1)*(fcut/2)).^2/(-2*log(GNyq(ii))));
        H = fspecial('gaussian', [N,1], alpha(2))*fspecial('gaussian', [1,N], alpha(1));
        Hd = H./max(H(:));
        h = fwind1(Hd,kaiser(N));
        I_MS_LR(:,:,ii) = imfilter(I_MS(:,:,ii),real(h),'replicate');
    end
    
    %%% Decimation MS
    if flag_decimate==1
        I_MS_LR = imresize(I_MS_LR,round([L1,L2]./ratio),'nearest');
    end

%% Test for non integer ratio (Incomplete and Unworking)
%     tol=10^(-5);
%     ratio_num=ratio; ratio_den=[1,1];
%     while ratio_den(1)<20 && rem(ratio_num(1),1)>tol
%         ratio_den(1)=ratio_den(1)+1;
%         ratio_num(1)=ratio(1)*ratio_den(1);
%     end
%     while ratio_den(2)<20 && rem(ratio_num(1),1)>tol
%         ratio_den(2)=ratio_den(2)+1;
%         ratio_num(2)=ratio(2)*ratio_den(1);
%     end
%     ratio_num=round(ratio_num);
%     if ratio_den==20, error('Ratio is not fractional'); end
%     
%     N = 41;
%     M = N*ratio_den+rem(ratio_den+1,2);
%     I_MS_LR = zeros([L1,L2,Nb]);
%     fcut = 1 ./ ratio_num;
%     
%     for ii = 1 : Nb
%         alpha = sqrt(((M-1).*(fcut/2)).^2/(-2*log(GNyq(ii))));
%         H = fspecial('gaussian', [M(2),1], alpha(2))*fspecial('gaussian', [1,M(1)], alpha(1));
%         Hd = H./max(H(:));
%         h = fwind1(Hd,kaiser(N));
%         
%         I_MS_LR(:,:,ii) = imfilter(I_MS(:,:,ii),real(h),'replicate');
%     end
    
end  

end


