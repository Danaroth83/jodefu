%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           GSA fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Gram-Schmidt Adaptive (GSA) algorithm.
% 
% Interface:
%           I_Fus_GSA = GSA(I_MS,I_PAN,I_MS_LR,ratio)
%
% Inputs:
%           I_MS:       MS image upsampled at PAN scale;
%           I_PAN:      PAN image;
%           I_MS_LR:    MS image;
%           ratio:      Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_GSA:  GSA pasharpened image.
% 
% References:
%           [Aiazzi07]   B. Aiazzi, S. Baronti, and M. Selva, “Improving component substitution Pansharpening through multivariate regression of MS+Pan
%                        data,” IEEE Transactions on Geoscience and Remote Sensing, vol. 45, no. 10, pp. 3230–3239, October 2007.
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_Fus_GSA,alpha,g] = GSA_local_HPM(I_MS,I_PAN,~,ratio,block_win,block_shift,flag_downsample)

if nargin<=4 || isempty(block_win), block_win=33; end
if nargin<=5 || isempty(block_shift), block_shift=1; end
if nargin<=6, flag_downsample=0; end

imageLR = double(I_MS);
imageHR = double(I_PAN);

%%% Remove means from imageLR
imageLR0 = zeros(size(I_MS));
for ii = 1 : size(I_MS,3), imageLR0(:,:,ii) = imageLR(:,:,ii) - mean2(imageLR(:,:,ii)); end

%%% Sintetic intensity
imageHR0 = imageHR - mean2(imageHR);
if flag_downsample==1
    imageHR0_LP = LowPass_Bspline(imageHR0,ratio);
else
    imageHR0_LP = LPfilterWithoutDec(imageHR0,ratio);
end
[~,alpha] = estimation_alpha(cat(3,imageLR0,ones(size(I_MS,1),size(I_MS,2))),imageHR0_LP,'local',block_win,block_shift);
alpha=imresize(alpha,block_shift);
alpha=padarray(alpha,[size(imageLR,1)-size(alpha,1),size(imageLR,2)-size(alpha,2)],'replicate','post');
I=sum(cat(3,imageLR0,ones([size(imageLR0,1),size(imageLR0,2)])).*alpha,3);
I0=I-mean2(I);
g=imageLR0./repmat(I0,[1,1,size(imageLR0,3)]);
delta=repmat(imageHR0-I0,[1,1,size(imageLR0,3)]);
I_Fus_GSA=imageLR0+g.*delta;

% Final Mean Equalization
for ii = 1 : size(I_MS,3)
    h = I_Fus_GSA(:,:,ii);
    I_Fus_GSA(:,:,ii) = h - mean2(h) + mean2(imageLR(:,:,ii));
end

end