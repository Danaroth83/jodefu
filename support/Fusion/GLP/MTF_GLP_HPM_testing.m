%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           MTF_GLP_HPM fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Modulation Transfer Function - Generalized Laplacian Pyramid (MTF-GLP) with High Pass Modulation (HPM) injection model algorithm. 
% 
% Interface:
%           I_Fus_MTF_GLP_HPM = MTF_GLP_HPM(I_PAN,I_MS,GNyq,ratio)
%
% Inputs:
%           I_PAN:              PAN image;
%           I_MS:               MS image upsampled at PAN scale;
%           GNyq:               MTF gain of MS sensor at Nyquist frequency
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_MTF_GLP_HPM:  MTF_GLP_HPM pansharpened image.
% 
% References:
%           [Aiazzi03]  B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “An MTF-based spectral distortion minimizing model for Pan-sharpening
%                       of very high resolution multispectral images of urban areas,” in Proceedings of URBAN 2003: 2nd GRSS/ISPRS Joint Workshop on
%                       Remote Sensing and Data Fusion over Urban Areas, 2003, pp. 90–94.
%           [Aiazzi06]  B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,”
%                       Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591–596, May 2006.
%           [Vivone14a] G. Vivone, R. Restaino, M. Dalla Mura, G. Licciardi, and J. Chanussot, “Contrast and error-based fusion schemes for multispectral
%                       image pansharpening,” IEEE Geoscience and Remote Sensing Letters, vol. 11, no. 5, pp. 930–934, May 2014.
%           [Vivone14b] G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                       IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_Fus_MTF_GLP_HPM,g,g_opt] = MTF_GLP_HPM_testing(I_PAN,I_MS,GNyq,ratio,I_GT)
% keyboard
imageHR = double(I_PAN);
I_MS = double(I_MS);

% %%% Equalization
imageHR = repmat(imageHR,[1 1 size(I_MS,3)]);

for ii = 1 : size(I_MS,3)  
 imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));  
end

% switch sensor
%     case 'QB' 
%         GNyq = [0.34 0.32 0.30 0.22]; % Band Order: B,G,R,NIR
%     case 'IKONOS'
%         GNyq = [0.26,0.28,0.29,0.28]; % Band Order: B,G,R,NIR
%     case 'GeoEye1'
%         GNyq = [0.23,0.23,0.23,0.23]; % Band Order: B,G,R,NIR
%     case 'WV2'
%         GNyq = [0.35 .* ones(1,7), 0.27];
%     case 'none'
%         if strcmp(tag,'WV2')
%             GNyq = 0.15 .* ones(1,8);
%         else
%             GNyq = 0.29 .* ones(1,size(I_MS,3)); %toolbox
% %             GNyq = 0.3 .* ones(1,size(I_MS,3));
%         end
% end

%%% MTF
N = 41;
PAN_LP = zeros(size(I_MS));
nBands = size(I_MS,3);
fcut = 1/ratio;
PSF_G = zeros(N,N,nBands);

for ii = 1 : nBands
    alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq(ii))));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    PSF_G(:,:,ii) = real(h);
    PAN_LP(:,:,ii) = imfilter(imageHR(:,:,ii),real(h),'replicate');
    t = imresize(PAN_LP(:,:,ii),1/ratio,'nearest');
    if (2^round(log2(ratio)) ~= ratio)
        PAN_LP(:,:,ii) = imresize(t,ratio);
    else
        PAN_LP(:,:,ii) = interp23tap(t,ratio);
    end
end

PAN_LP = double(PAN_LP);
g=(I_MS ./ (PAN_LP + eps));
g_opt=(I_GT ./ (imageHR + eps));
I_Fus_MTF_GLP_HPM = I_MS .* (imageHR ./ (PAN_LP + eps));

end