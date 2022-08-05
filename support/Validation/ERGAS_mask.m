%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Erreur Relative Globale Adimensionnelle de Synthèse (ERGAS).
% 
% Interface:
%           ERGAS_index = ERGAS(I1,I2,ratio)
%
% Inputs:
%           I1:             First multispectral image;
%           I2:             Second multispectral image;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Integer value.
% 
% Outputs:
%           ERGAS_index:    ERGAS index.
% References:
%           [Ranchin00]     T. Ranchin and L. Wald, “Fusion of high spatial and spectral resolution images: the ARSIS concept and its implementation,”
%                           Photogrammetric Engineering and Remote Sensing, vol. 66, no. 1, pp. 49–61, January 2000.
%           [Vivone14]      G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                           IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ERGAS_index,ERGAS_map] = ERGAS_mask(I1,I2,ratio,mask)

if nargin<=3, mask=ones(size(I1,1),size(I1,2)); end
if nargout>1, ERGAS_map=zeros(size(I1,1),size(I1,2)); end

Nc=max(mask(:));
Nb=size(I1,3);
ERGAS_index=zeros(1,Nc);

for nc=1:Nc
    mask_temp=(mask==nc);
    Nm=nnz(mask_temp);

    I1_p=zeros(Nm,Nb);
    I2_p=zeros(Nm,Nb);
    for iLR=1:Nb
        I1_v = I1(:,:,iLR);
        I2_v = I2(:,:,iLR);
        I1_p(:,iLR) = double(I1_v(mask_temp));
        I2_p(:,iLR) = double(I2_v(mask_temp));
    end

    Err=I1_p-I2_p;
    for iLR=1:Nb
        ERGAS_index(nc) = ERGAS_index(nc)+mean(Err(:,iLR).^2)/(mean(I1_p(:,iLR)))^2;   
    end
    ERGAS_index(nc) = (100/ratio) * sqrt((1/Nb) * ERGAS_index(nc));

    if nargout>1
        ERGAS_map_p=zeros(Nm,1);
        for iLR=1:Nb
            ERGAS_map_p = ERGAS_map_p+(Err(:,iLR).^2)/(mean(I1_p(:,iLR)))^2;   
        end
        ERGAS_map_v = (100/ratio) * sqrt((1/Nb) * ERGAS_map_p);
        ERGAS_map(mask_temp) = ERGAS_map_v;
    end
end

ERGAS_map(mask<=0)=NaN;

end