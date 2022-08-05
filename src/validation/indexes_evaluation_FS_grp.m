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
function [D_lambda,D_S,QNR_index,SAM_index,sCC,D_S2,SAM_index2,sCC_Sobel] = indexes_evaluation_FS_grp(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,band_overlap_cell,LR_method)
if nargin<=11
    LR_method='MTF';
end
switch LR_method
    case 1
        LR_method='resize';
    case 0
        LR_method='MTF';
end

if th_values
    I_F(I_F > 2^L) = 2^L;
    I_F(I_F < 0) = 0;
end

QNR_index=zeros(1,size(I_PAN,3));
D_lambda=zeros(1,size(I_PAN,3));
D_S=zeros(1,size(I_PAN,3));
D_S2=zeros(1,size(I_PAN,3));
for ii=1:size(I_PAN,3)
    % [QNR_index(ii),D_lambda(ii),D_S(ii)]= QNR(I_F,I_MS,I_PAN,sensor,ratio,Qblocks_size);
    if nargout>=6
        [QNR_index(ii),D_lambda(ii),D_S(ii),D_S2(ii)]= QNR_mod(I_F(:,:,band_overlap_cell{ii}),I_MS_LR(:,:,band_overlap_cell{ii}),I_MS(:,:,band_overlap_cell{ii}),I_PAN(:,:,ii),GNyq_PAN(ii),ratio,Qblocks_size,LR_method);
        D_S2(ii) = D_S2(ii) * length(band_overlap_cell{ii});
    else
        [QNR_index(ii),D_lambda(ii),D_S(ii)]= QNR_mod(I_F(:,:,band_overlap_cell{ii}),I_MS_LR(:,:,band_overlap_cell{ii}),I_MS(:,:,band_overlap_cell{ii}),I_PAN(:,:,ii),GNyq_PAN(ii),ratio,Qblocks_size,LR_method);
    end
    QNR_index(ii) = QNR_index(ii) * length(band_overlap_cell{ii});
    D_lambda(ii) = D_lambda(ii) * length(band_overlap_cell{ii});
    D_S(ii) = D_S(ii) * length(band_overlap_cell{ii});
end
QNR_index=sum(QNR_index)/length(cell2mat(band_overlap_cell));
D_lambda=sum(D_lambda)/length(cell2mat(band_overlap_cell));
D_S=sum(D_S)/length(cell2mat(band_overlap_cell));
if nargout>=6
    D_S2=sum(D_S2)/length(cell2mat(band_overlap_cell));
end

% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

SAM_index = SAM(I_MS_LR,I_F_D);
if nargout>=7
    SAM_index2 = SAM(I_MS,I_F);
end

sCC=zeros(1,size(I_PAN,3));
sCC_Sobel=zeros(1,size(I_PAN,3));
for ii=1:size(I_PAN,3)
    I_PAN_temp = repmat(I_PAN(:,:,ii),[1 1 length(band_overlap_cell{ii})]);

    Im_Lap_PAN = zeros(size(I_PAN_temp));
    for idim=1:size(I_PAN_temp,3)
        % Im_Lap_PAN(:,:,idim)= imfilter(I_PAN_temp(:,:,idim),fspecial('sobel'));
        Im_Lap_PAN(:,:,idim)= imfilter(I_PAN_temp(:,:,idim),fspecial('laplacian'));
    end

    % I_F(I_F > 2^L) = 2^L;
    % I_F(I_F < 0) = 0;
    
    I_F_temp=I_F(:,:,band_overlap_cell{ii});
    Im_Lap_F = zeros(size(I_F_temp));
    for idim=1:length(band_overlap_cell{ii})
        % Im_Lap_F(:,:,idim)= imfilter(I_F_temp(:,:,idim),fspecial('sobel'));
        Im_Lap_F(:,:,idim)= imfilter(I_F_temp(:,:,idim),fspecial('laplacian'));
    end
    % keyboard
    %  Im_Lap_F(Im_Lap_F > 2^16) = 2^16;
    %  Im_Lap_F(Im_Lap_F <- 2^16) = -2^16;

    sCC(ii)=sum(sum(sum(Im_Lap_PAN.*Im_Lap_F)));
    sCC(ii)=sCC(ii)/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))));
    sCC(ii)=sCC(ii)/sqrt(sum(sum(sum((Im_Lap_F.^2)))));
    
    sCC(ii)=sCC(ii)*length(band_overlap_cell{ii});
    
    if nargout>=8
        for idim=1:size(I_PAN_temp,3)
            Im_Lap_PAN(:,:,idim)= imfilter(I_PAN_temp(:,:,idim),fspecial('sobel'));
        end
        for idim=1:length(band_overlap_cell{ii})
            Im_Lap_F(:,:,idim)= imfilter(I_F_temp(:,:,idim),fspecial('sobel'));
        end
        sCC_Sobel(ii)=sum(sum(sum(Im_Lap_PAN.*Im_Lap_F)));
        sCC_Sobel(ii)=sCC_Sobel(ii)/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))));
        sCC_Sobel(ii)=sCC_Sobel(ii)/sqrt(sum(sum(sum((Im_Lap_F.^2)))));

        sCC_Sobel(ii)=sCC_Sobel(ii)*length(band_overlap_cell{ii});
    end
end
sCC=sum(sCC)/length(cell2mat(band_overlap_cell));
sCC_Sobel=sum(sCC_Sobel)/length(cell2mat(band_overlap_cell));

end