%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Spectral Angle Mapper (SAM).
% 
% Interface:
%           [SAM_index,SAM_map] = SAM(I1,I2)
%
% Inputs:
%           I1:         First multispectral image;
%           I2:         Second multispectral image.
% 
% Outputs:
%           SAM_index:  SAM index;
%           SAM_map:    Image of SAM values.
% 
% References:
%           [Yuhas92]   R. H. Yuhas, A. F. H. Goetz, and J. W. Boardman, "Discrimination among semi-arid landscape endmembers using the Spectral Angle Mapper (SAM) algorithm," 
%                       in Proceeding Summaries 3rd Annual JPL Airborne Geoscience Workshop, 1992, pp. 147–149.
%           [Vivone14]  G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                       IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SAM_index,SAM_map] = SAM_mask(I1,I2,mask)

if nargin<=2, mask=ones(size(I1,1),size(I1,2)); end

Nc=max(mask(:));
SAM_index=zeros(1,Nc);
SAM_map = zeros(size(I1,1),size(I1,2));
Nb=size(I1,3);

for nc=1:Nc
    mask_temp=(mask==nc);
    Nm=nnz(mask_temp);

    I1_masked=zeros(Nm,Nb);
    I2_masked=zeros(Nm,Nb);

    for ii=1:Nb
        I1_temp=I1(:,:,ii);
        I2_temp=I2(:,:,ii);
        I1_masked(:,ii)=I1_temp(mask_temp);
        I2_masked(:,ii)=I2_temp(mask_temp);
    end

    prod_scal = dot(I1_masked,I2_masked,2); 
    norm_orig = dot(I1_masked,I1_masked,2);
    norm_fusa = dot(I2_masked,I2_masked,2);
    prod_norm = sqrt(norm_orig.*norm_fusa);
    prod_map = prod_norm;
    prod_map(prod_map==0)=eps;
    SAM_map(mask_temp) = real(acos(prod_scal./prod_map))*180/pi;
    z=find(prod_norm==0);
    prod_scal(z)=[];prod_norm(z)=[];
    angolo = sum(sum(acos(prod_scal./prod_norm)))/(size(prod_norm,1));
    SAM_index(nc) = real(angolo)*180/pi;
end
SAM_map(mask<=0)=NaN;

end