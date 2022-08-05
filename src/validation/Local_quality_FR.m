function [SAM_index,SAM_map, sCC,sCC_map, Im_Lap_F,equivPAN,Q_PAN,Q_PAN_map,Im_Edges]=Local_quality_FR(I_F,I_MS_LR,I_MS,I_PAN,Q_blocks_size,Q_shift,ratio,flag_cut_bounds,dim_cut,thvalues,sensor,im_tag,L)

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
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           tag:                Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number;
%
% Outputs:
%           Q_index, Q_map:         Q index and map; 
%           SAM_index, SAM_map:     Spectral Angle Mapper (SAM) index and map;
%           ERGAS_index, ERGAS_map: Erreur Relative Globale Adimensionnelle de Synthèse (ERGAS) index and map;
%           Q2n_index, Q2n_map:     Q2n indexand map.
%           The maps have size [N1-(2*dim_cut-1)) x N2-(2*dim_cut-1)]       
%           The Q_map has size [N1-(2*dim_cut-1)-Q_blocks_size/2) x N2-(2*dim_cut-1)-Q_blocks_size/2)] 
%                                   since it is calculated only on valid values (without convolution aliasing)

% keyboard
% 
% if flag_cut_bounds
%     I_MS = I_MS(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
%     I_PAN = I_PAN(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
%     I_F = I_F(dim_cut:end-dim_cut,dim_cut:end-dim_cut,:);
%     for iband=1:size(I_MS,3)
%         I_MS_LR(:,:,iband) = imresize(I_MS_LR(:,:,iband),[size(I_MS,1)*ratio,size(I_MS,2)*ratio]);
%     end
% end

if thvalues
    I_F(I_F > 2^L) = 2^L;
    I_F(I_F < 0) = 0;
end


I_F = double(I_F);
I_PAN = double(I_PAN);
I_MS_LR = double(I_MS_LR);

cd Quality_Indices

%%% Reduced scale version
I_F_LP = MTF(I_F,sensor,im_tag,ratio);
I_F_HP=I_F-I_F_LP;
I_F_D = imresize(I_F_LP,1/ratio,'nearest');
%%%Estimate at full scale 
%%Direct
PAN_LP_estim = I_PAN;
%%Con Imresize
% PAN_LP_estim=imresize(I_PAN,1/ratio);  
% if  (2^round(log2(ratio)) ~= ratio)
%    PAN_LP_estim = imresize(PAN_LP_estim,[size(I_PAN,1),size(I_PAN,1)]);
% else
%     PAN_LP_estim = interp23tapGeneral(PAN_LP_estim,ratio);
% end
% 
% Old
% alpha(1,1,:) = estimation_alpha(cat(3,I_F,ones(size(I_F,1),size(I_F,2))),PAN_LP_estim ,'global');
% equivPAN=sum(cat(3,I_F,ones(size(I_F,1),size(I_F,2))) .* repmat(alpha,[size(I_F,1) size(I_F,2) 1]),3);


%%% Estimate at reduced scale
%%Con Imresize
% PAN_LP_estim=imresize(I_PAN,1/ratio);  
%%Con MTF_PAN
% PAN_LP_estim1=MTF_PAN(I_PAN,sensor,ratio);
%  PAN_LP_estim=imresize(PAN_LP_estim1,1/ratio,'nearest');  
% alpha(1,1,:) = estimation_alpha(cat(3,I_F_D,ones(size(I_F_D,1),size(I_F_D,2))),PAN_LP_estim,'global');
I_PAN_LR=imresize(I_PAN,1/ratio);
I_MS_LR_C = reshape(I_MS_LR,[size(I_MS_LR,1)*size(I_MS_LR,2),size(I_MS_LR,3)]);
I_PAN_LR_C = I_PAN_LR(:);
alpha=pinv([I_MS_LR_C,ones(size(I_MS_LR_C,1),1)])*I_PAN_LR_C;

%Reconstruction
I_F_C = reshape(I_F,[size(I_F,1)*size(I_F,2),size(I_F,3)]);
I_F_PAN_C=[I_F_C,ones(size(I_F_C,1),1)]*alpha;
equivPAN = reshape(I_F_PAN_C,[size(I_PAN,1),size(I_PAN,2)]);

% keyboard
[Q_PAN,Q_PAN_map]=img_qi(I_PAN,equivPAN,Q_blocks_size);



imageHR = repmat(I_PAN,[1 1 size(I_MS_LR,3)]);

Im_Lap_PAN = zeros(size(imageHR));
for idim=1:size(imageHR,3),
    %     Im_Lap_PAN(:,:,idim)= imfilter(imageHR(:,:,idim),fspecial('sobel'));
    Im_Lap_PAN_y(:,:,idim)= imfilter(imageHR(:,:,idim),fspecial('sobel'));
    Im_Lap_PAN_x(:,:,idim)= imfilter(imageHR(:,:,idim),fspecial('sobel')');
    Im_Lap_PAN(:,:,idim) = sqrt(Im_Lap_PAN_y(:,:,idim).^2+Im_Lap_PAN_x(:,:,idim).^2);
end

Im_Lap_F = zeros(size(I_F));
for idim=1:size(I_MS_LR,3),
    Im_Lap_F1(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
    Im_Lap_F2(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel')');
    Im_Lap_F(:,:,idim)=Im_Lap_F1(:,:,idim).^2+Im_Lap_F2(:,:,idim).^2;
end

Im_Edges = zeros(size(I_F));
for idim=1:size(I_MS_LR,3),
    % Im_Lap_Edges(:,:,idim)= edge(I_F_HP(:,:,idim),'canny');
    Im_Edges(:,:,idim)= edge(I_F(:,:,idim),'canny');
end
sCC=sum(sum(sum(Im_Lap_PAN.*Im_Lap_F)));
sCC=sCC/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))));
sCC=sCC/sqrt(sum(sum(sum((Im_Lap_F.^2)))));
% sCC_map=sum(Im_Lap_PAN.*Im_Lap_F,3);
sCC_map=Im_Lap_PAN.*Im_Lap_F;

sCC_map=sCC_map/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))))/sqrt(sum(sum(sum((Im_Lap_F.^2)))));


% keyboard
[SAM_index,SAM_map] = SAM(I_MS_LR,I_F_D);

cd ..

% 
% for idim=1:size(I_F,3),
%     I_GT_band=I_GT(:,:,idim);
%     I_F_band=I_F(:,:,idim);
%     [Q_index_band(idim),Q_map_band(:,:,idim)] = img_qi(I_GT_band,I_F_band, Q_blocks_size);
%     ERGAS_map_band(:,:,idim) = (I_GT_band-I_F_band).^2/(mean2(I_GT_band))^2;
%     ERGAS_band(idim)=mean2((I_GT_band-I_F_band).^2)/(mean2(I_GT_band))^2;
% end
% % keyboard
% ERGAS_map = mean(ERGAS_map_band,3);
% ERGAS = (100/ratio) * sqrt((1/size(I_F,3)) * sum(ERGAS_band));
% 
% [SAM_index,SAM_map] = SAM(I_GT,I_F);
% [Q2n_index, Q2n_map] = q2n(I_GT,I_F, Q_blocks_size, Q_shift);
% % keyboard
% Q_index=mean(Q_index_band);
% Q_map=mean(Q_map_band,3);
% cd ..
