%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%           MTF filters the image I_MS using a Gaussin filter matched with the Modulation Transfer Function (MTF) of the MultiSpectral (MS) sensor.
%
% Interface:
%           I_Filtered = MTF(I_MS,sensor,tag,ratio)
%
% Inputs:
%           I_MS:           MS image;
%           sensor:         String for type of sensor (e.g. 'WV2', 'IKONOS');
%           tag:            Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number
%                           in the latter case.
%           ratio:          Scale ratio between MS and PAN.
%
% Outputs:
%           I_Filtered:     Output filtered MS image.
%
% References:
%           [Aiazzi06]   B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, �MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,�
%                        Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591�596, May 2006.
%           [Lee10]      J. Lee and C. Lee, �Fast and efficient panchromatic sharpening,� IEEE Transactions on Geoscience and Remote Sensing, vol. 48, no. 1,
%                        pp. 155�163, January 2010.
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, �A Critical Comparison Among Pansharpening Algorithms�,
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_Filtered = MTF(I_MS,sensor,tag,ratio)

switch sensor
    case 'QB'
        GNyq = [0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
    case 'IKONOS'
        GNyq = [0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
    case 'GeoEye1'
        GNyq = [0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
    case 'WV2'
        GNyq = [0.35 .* ones(1,7), 0.27];
    case {'WV3','WV3_4bands'}
        GNyq = [0.325 0.355 0.360 0.350 0.365 0.360 0.335 0.315];
        if strcmp(sensor,'WV3_4bands')
            GNyq = GNyq([2 3 5 7]);
        end
    case {'HYP','HYP_14_33','HYP_16_31'}
        flag_resize_new = 2; % MTF usage
        %VNIR
        MTF_MS(1:21)=0.27;
        MTF_MS(22:41)=0.28;
        MTF_MS(42:49)=0.26;
        MTF_MS(50:70)=0.26;
        %SWIR
        MTF_MS(71:100)=0.30;
        MTF_MS(101:130)=0.30;
        MTF_MS(131:177)=0.27;
        MTF_MS(177:242)=0.27;
        if strcmp(sensor,'HYP_14_33')
            GNyq = MTF_MS(14:33);
        elseif strcmp(sensor,'HYP_16_31')
            GNyq = MTF_MS(16:31);
        else
            GNyq = MTF_MS;
        end
    case 'none'
        GNyq = 0.29 .* ones(1,size(I_MS,3));
end


%%% MTF

N = 41;
I_MS_LP = zeros(size(I_MS));
nBands = size(I_MS,3);
fcut = 1/ratio;

for ii = 1 : nBands
    alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq(ii))));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    I_MS_LP(:,:,ii) = imfilter(I_MS(:,:,ii),real(h),'replicate');
end

I_Filtered= double(I_MS_LP);
end