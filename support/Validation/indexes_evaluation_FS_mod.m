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
function [D_lambda,D_S,QNR_index,SAM_index,sCC,D_S2,SAM_index2,sCC_Sobel] = indexes_evaluation_FS_mod(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,LR_method)
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

if nargout<=5
    % [QNR_index,D_lambda,D_S]= QNR(I_F,I_MS,I_PAN,sensor,ratio,Qblocks_size);
    [QNR_index,D_lambda,D_S]= QNR_mod(I_F,I_MS_LR,I_MS,I_PAN,GNyq_PAN,ratio,Qblocks_size,LR_method);
else
    [QNR_index,D_lambda,D_S,D_S2]= QNR_mod(I_F,I_MS_LR,I_MS,I_PAN,GNyq_PAN,ratio,Qblocks_size,LR_method);
end

%     I_F(I_F > 2^L) = 2^L;
%     I_F(I_F < 0) = 0;

I_PAN = repmat(I_PAN,[1 1 size(I_MS,3)]);

% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

SAM_index = SAM(I_MS_LR,I_F_D); %toolbox
if nargout>=7
    SAM_index2 = SAM(I_MS,I_F); % HS toolbox
end

Im_Lap_PAN = zeros(size(I_PAN));
for idim=1:size(I_PAN,3)
    Im_Lap_PAN(:,:,idim)= imfilter(I_PAN(:,:,idim),fspecial('laplacian'));
end

Im_Lap_F = zeros(size(I_F));
for idim=1:size(I_MS_LR,3)
    Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('laplacian'));
end

% keyboard
%  Im_Lap_F(Im_Lap_F > 2^16) = 2^16;
%  Im_Lap_F(Im_Lap_F <- 2^16) = -2^16;

sCC=sum(sum(sum(Im_Lap_PAN.*Im_Lap_F)));
sCC=sCC/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))));
sCC=sCC/sqrt(sum(sum(sum((Im_Lap_F.^2)))));

if nargout>=8
    for idim=1:size(I_PAN,3)
        Im_Lap_PAN(:,:,idim)= imfilter(I_PAN(:,:,idim),fspecial('sobel'));
    end
    for idim=1:size(I_MS_LR,3)
        Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
    end
    sCC_Sobel=sum(sum(sum(Im_Lap_PAN.*Im_Lap_F)));
    sCC_Sobel=sCC_Sobel/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))));
    sCC_Sobel=sCC_Sobel/sqrt(sum(sum(sum((Im_Lap_F.^2)))));
end

end