%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           ATWT fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the A Trous Wavelet Transform (ATWT) and additive injection model algorithm.
% 
% Interface:
%           I_Fus_ATWT = ATWT(I_MS,I_PAN,ratio)
%
% Inputs:
%           I_MS:       MS image upsampled at PAN scale;
%           I_PAN:      PAN image;
%           ratio:      Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_ATWT: ATWT pasharpened image.
% 
% References:
%           [Nunez99]       J. Nunez, X. Otazu, O. Fors, A. Prades, V. Pala, and R. Arbiol, “Multiresolution-based image fusion with additive wavelet
%                           decomposition,” IEEE Transactions on Geoscience and Remote Sensing, vol. 37, no. 3, pp. 1204–1211, May 1999.
%           [Vivone14a]     G. Vivone, R. Restaino, M. Dalla Mura, G. Licciardi, and J. Chanussot, “Contrast and error-based fusion schemes for multispectral
%                           image pansharpening,” IEEE Geoscience and Remote Sensing Letters, vol. 11, no. 5, pp. 930–934, May 2014.
%           [Vivone14b]     G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                           IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_Fus_ATWT = ATWT_mod(I_MS,I_PAN,ratio,I_PAN_eq)

if nargin<=3
    %%% Equalization for each band
    mean2_I_PAN=mean2(I_PAN);
    std2_I_PAN=std2(I_PAN);
    I_PAN = repmat(I_PAN,[1 1 size(I_MS,3)]);
    for ii = 1 : size(I_MS,3)    
      I_PAN(:,:,ii) = (I_PAN(:,:,ii) - mean2_I_PAN).*(std2(I_MS(:,:,ii))./std2_I_PAN) + mean2(I_MS(:,:,ii));  
    end
else
    I_PAN=I_PAN_eq;
end

[Height,Width,Bands]=size(I_MS);
I_Fus_ATWT=zeros(Height,Width,Bands,'double');


%%% Different w.r.t. [Nunez99] in which the equalization is done
%%% w.r.t. the Luminance of the MS image



% h=[1 4 6 4 1 ]/16;
% g=[0 0 1 0 0 ]-h;
% htilde=[ 1 4 6 4 1]/16;
% gtilde=[ 0 0 1 0 0 ]+htilde;
% h=sqrt(2)*h;
% g=sqrt(2)*g;
% htilde=sqrt(2)*htilde;
% gtilde=sqrt(2)*gtilde;
% WF={h,g,htilde,gtilde};

for i=1:Bands    
    WT = swt2_Bspline(I_PAN(:,:,i),ratio);
    
    
    StepDetails = I_PAN(:,:,i) - iswt2_Bspline(WT,ratio);
    
%%%%%%%%% OLD (as in the article [Nunez99]). Lower performances.
%     sINI = WT.sizeINI;
%     
%     StepDetails = zeros(sINI);
%     
%     for ii = 2 : numel(WT.dec)    
%         h = WT.dec{ii};
%         h = imcrop(h,[(size(h,1) - sINI(1))/2 + 1,(size(h,2) - sINI(2))/2 + 1, sINI(1) - 1, sINI(2) - 1]);
%         StepDetails = StepDetails + h; 
%     end
%%%%%%%%%%%%%%%%%%%
    
    I_Fus_ATWT(:,:,i) = StepDetails + I_MS(:,:,i);
end

end