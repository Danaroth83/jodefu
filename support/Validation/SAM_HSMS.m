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
function [SAM_index,SAM_map] = SAM_HSMS(I1,I2MS)


[M,N,Nb] = size(I1);
Nb_MS = size(I2MS,3);
I1dotI1=I1.*I1;
norm_I1 = sum(I1dotI1,3);
norm_I1_rep=repmat(norm_I1,[1,1,Nb_MS]);
SAM_map=zeros(M,N,Nb_MS,Nb);
SAM_index=zeros(Nb_MS,Nb);
for ii=1:Nb
    I2MS_temp=I2MS(:,:,:,ii);
    norm_I1ii_rep=repmat(I1dotI1(:,:,ii),[1,1,Nb_MS]);
    norm_I2_rep=norm_I1_rep-norm_I1ii_rep+I2MS_temp.*I2MS_temp;
    norm_I1I2_rep=norm_I1_rep-norm_I1ii_rep+repmat(I1(:,:,ii),[1,1,Nb_MS]).*I2MS_temp;
    prod_norm=sqrt(norm_I2_rep.*norm_I1_rep);
    prod_norm(prod_norm==0)=eps;
    SAM_temp=acos(norm_I1I2_rep./prod_norm);
    SAM_map(:,:,:,ii)=SAM_temp;
    SAM_index(:,ii)=real(sum(sum(SAM_temp,2),1))/(N*M)*180/pi;
end

end