function [Q2n_index, Q2n_map,Q_index,Q_map,Q_map_band,ERGAS,ERGAS_map,...
    sCC_index,sCC_map,SAM_index,SAM_map]=Local_quality(I_F,I_GT,Q_blocks_size,Q_shift,ratio,flag_cut_bounds,dim_cut,thvalues)

% Inputs:
%           I_F:                    Fused Image [N1 x N2 x N];
%           I_GT:                   Ground-Truth image [N1 x N2 x N];
%           ratio:                  Scale ratio between MS and PAN. Pre-condition: Integer value;
%           L:                      Image radiometric resolution; 
%           Q_blocks_size:          Block size of the Q-index locally applied;
%           Q_shift:                Discretization of the Q map of the Q-index locally applied;
%           flag_cut_bounds:        Cut the boundaries of the viewed Panchromatic image;
%           dim_cut:                Define the dimension of the boundary cut;
%           th_values:              Flag. If th_values == 1, apply an hard threshold to the dynamic range.
%
% Outputs:
%           Q_index, Q_map:         Q index and map; 
%           SAM_index, SAM_map:     Spectral Angle Mapper (SAM) index and map;
%           sCC_index, sCC_map:     spatial Correlation Coefficient (sCC) index and map;
%           ERGAS_index, ERGAS_map: Erreur Relative Globale Adimensionnelle de Synthèse (ERGAS) index and map;
%           Q2n_index, Q2n_map:     Q2n indexand map.
%           The maps have size [N1-(2*dim_cut-1)) x N2-(2*dim_cut-1)]       
%           The Q_map has size [N1-(2*dim_cut-1)-Q_blocks_size/2) x N2-(2*dim_cut-1)-Q_blocks_size/2)] 
%                                   since it is calculated only on valid values (without convolution aliasing)

if flag_cut_bounds
    I_GT = I_GT(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
    I_F = I_F(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
end

if thvalues
    I_F(I_F > 2^L) = 2^L;
    I_F(I_F < 0) = 0;
end

I_F = double(I_F);
I_GT = double(I_GT);

cd Quality_Indices
for idim=1:size(I_F,3),
    I_GT_band=I_GT(:,:,idim);
    I_F_band=I_F(:,:,idim);
    [Q_index_band(idim),Q_map_band(:,:,idim)] = img_qi(I_GT_band,I_F_band, Q_blocks_size);
    ERGAS_map_band(:,:,idim) = (I_GT_band-I_F_band).^2/(mean2(I_GT_band))^2;
    ERGAS_band(idim)=mean2((I_GT_band-I_F_band).^2)/(mean2(I_GT_band))^2;
end
% keyboard
ERGAS_map = mean(ERGAS_map_band,3);
ERGAS = (100/ratio) * sqrt((1/size(I_F,3)) * sum(ERGAS_band));

[SAM_index,SAM_map] = SAM(I_GT,I_F);
[Q2n_index, Q2n_map] = q2n(I_GT,I_F, Q_blocks_size, Q_shift);
% keyboard
Q_index=mean(Q_index_band);
Q_map=mean(Q_map_band,3);

%%% sCC
% Im_Lap_F = zeros(size(I_F));
% for idim=1:size(I_F,3),
%     Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
% end
% Im_Lap_GT = zeros(size(I_GT));
% for idim=1:size(I_GT,3),
%     Im_Lap_GT(:,:,idim)= imfilter(I_GT(:,:,idim),fspecial('sobel'));
% end

Im_Lap_GT = zeros(size(I_GT));
for idim=1:size(I_GT,3),
    %     Im_Lap_PAN(:,:,idim)= imfilter(imageHR(:,:,idim),fspecial('sobel'));
    Im_Lap_GT_y(:,:,idim)= imfilter(I_GT(:,:,idim),fspecial('sobel'));
    Im_Lap_GT_x(:,:,idim)= imfilter(I_GT(:,:,idim),fspecial('sobel')');
    Im_Lap_GT(:,:,idim) = sqrt(Im_Lap_GT_y(:,:,idim).^2+Im_Lap_GT_x(:,:,idim).^2);
end

Im_Lap_F = zeros(size(I_F));
for idim=1:size(I_F,3),
    Im_Lap_F1(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
    Im_Lap_F2(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel')');
    Im_Lap_F(:,:,idim)=Im_Lap_F1(:,:,idim).^2+Im_Lap_F2(:,:,idim).^2;
end

sCC = sum(sum(sum(Im_Lap_GT.*Im_Lap_F)));
sCC = sCC/sqrt(sum(sum(sum((Im_Lap_GT.^2)))));
sCC_index = sCC/sqrt(sum(sum(sum((Im_Lap_F.^2)))));
sCC_map=Im_Lap_GT.*Im_Lap_F;
sCC_map=sCC_map/sqrt(sum(sum(sum((Im_Lap_GT.^2)))))/sqrt(sum(sum(sum((Im_Lap_F.^2)))));


cd ..
