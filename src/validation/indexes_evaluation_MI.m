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
function [D_lambda,D_S,QNR_index,SAM_index,sCC,D_lambda_Zhou,D_S_Zhou,D_lambda_Khan,D_S_Khan] = indexes_evaluation_MI(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,LR_method)
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

[QNR_index,D_lambda,D_S]= QNR_mod(I_F,I_MS_LR,I_MS,I_PAN,GNyq_PAN,ratio,Qblocks_size,LR_method);

%     I_F(I_F > 2^L) = 2^L;
%     I_F(I_F < 0) = 0;

% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

SAM_index = SAM(I_MS_LR,I_F_D); %toolbox


Im_Lap_PAN= imfilter(I_PAN,fspecial('sobel'));
Im_Lap_PAN= repmat(Im_Lap_PAN,[1,1,size(I_MS,3)]);

Im_Lap_F = zeros(size(I_F));
for idim=1:size(I_MS_LR,3)
    Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
end

% keyboard
%  Im_Lap_F(Im_Lap_F > 2^16) = 2^16;
%  Im_Lap_F(Im_Lap_F <- 2^16) = -2^16;

sCC=SCC(Im_Lap_PAN,Im_Lap_F);

err=abs(I_F-I_MS);
D_lambda_Zhou=mean(err(:));
D_S_Zhou=1-sCC;

%Khan index
[D_lambda_Khan, D_S_Khan]=Khan_mod(I_F,I_MS,I_PAN,I_MS_LR,GNyq,ratio,Qblocks_size,LR_method,GNyq_PAN);

end