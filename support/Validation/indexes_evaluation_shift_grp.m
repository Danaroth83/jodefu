%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Full resolution quality indexes. 
% 
% Interface:
%           [D_lambda,D_S,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,sensor,tag,ratio)
%
% Inputs:
%           I_F:                Fused Image;
%           I_MS_LR:            MS image;
%           I_PAN:              Panchromatic image;
%           L:                  Image radiometric resolution; 
%           th_values:          Flag. If th_values == 1, apply an hard threshold to the dynamic range.
%           I_MS:               MS image upsampled to the PAN size;
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           tag:                Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           D_lambda:           D_lambda index;
%           D_S:                D_S index;
%           QNR_index:          QNR index;
%           SAM_index:          Spectral Angle Mapper (SAM) index between fused and MS image;
%           sCC:                spatial Correlation Coefficient between fused and PAN images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D_S,SAM_index,sCC,UIQI] = indexes_evaluation_shift_grp(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,band_overlap_cell,LR_method,S)
if nargin<=11 || isempty(S)
    S=1;
end
if nargin<=10 || isempty(LR_method)
    LR_method='resize';
end
switch LR_method
    case {1,'imresize'}
        LR_method='resize';
    case {0,'MTFfiltered'}
        LR_method='MTF';
end

if th_values
    I_F(I_F > 2^L) = 2^L;
    I_F(I_F < 0) = 0;
end

cd Quality_Indices_HS

% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);
SAM_index = SAM(I_MS_LR,I_F_D); %toolbox

D_S=0;
UIQI=0;
sCC=0;
for ii=1:size(I_PAN,3)
    band_overlap_temp=band_overlap_cell{ii};
    I_F_temp=I_F(:,:,band_overlap_temp);
    I_PAN_temp=I_PAN(:,:,ii);
    I_MS_temp=I_MS(:,:,band_overlap_temp);
    I_MS_LR_temp=I_MS_LR(:,:,band_overlap_temp);
    I_PAN_LR_F = LPF_filter(I_PAN_temp,ratio,LR_method,false,GNyq_PAN(ii));
    I_PAN_LR = LPF_filter(I_PAN_temp,ratio,LR_method,true,GNyq_PAN(ii));
    q=1;
    D_S_temp=spatial_distortion_new_mod(I_F_temp,I_PAN_temp,I_MS_temp,I_PAN_LR_F,q,Qblocks_size,S);
    
    Im_Lap_PAN=imfilter(I_PAN_temp,fspecial('sobel'));
    alpha=estimation_alpha(cat(3,ones(size(I_MS_LR_temp,1),size(I_MS_LR_temp,2)),I_MS_LR_temp),I_PAN_LR,'global');
    I_Sim=sum(cat(3,ones(size(I_F_temp,1),size(I_F_temp,2)),I_F_temp).*repmat(reshape(alpha,[1,1,length(alpha)]),[size(I_F_temp,1),size(I_F_temp,2),1]),3);
    Im_Lap_Sim=imfilter(I_Sim,fspecial('sobel'));
    sCC_temp=SCC(Im_Lap_Sim,Im_Lap_PAN,'local',Qblocks_size,S);
    UIQI_temp=img_qi_shift(Im_Lap_Sim,Im_Lap_PAN,Qblocks_size,S);
    
    D_S=D_S+D_S_temp*length(band_overlap_temp);
    sCC=sCC+sCC_temp*length(band_overlap_temp);
    UIQI=UIQI+UIQI_temp*length(band_overlap_temp);
end
D_S=D_S/length(cell2mat(band_overlap_cell));
sCC=sCC/length(cell2mat(band_overlap_cell));
UIQI=UIQI/length(cell2mat(band_overlap_cell));

cd ..

end