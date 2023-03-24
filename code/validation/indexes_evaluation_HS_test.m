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
function [D_lambda,D_S,QNR_index,SAM_index,sCC,D_S_Sob,D_S_HR_Sob,sCC_Sim_Sob,sCC_Sim_Lap,sCC_Sim_Sob_l,UIQI_Sim_Sob_l] = indexes_evaluation_HS_test(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,LR_method)
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

% [QNR_index,D_lambda,D_S]= QNR(I_F,I_MS,I_PAN,sensor,ratio,Qblocks_size);
[QNR_index,D_lambda,D_S]= QNR_mod(I_F,I_MS_LR,I_MS,I_PAN,GNyq_PAN,ratio,Qblocks_size,LR_method);

%     I_F(I_F > 2^L) = 2^L;
%     I_F(I_F < 0) = 0;


% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

SAM_index = SAM(I_MS_LR,I_F_D); %toolbox
% if nargout>=7
%     SAM_index2 = SAM(I_MS,I_F); % HS toolbox
% end


% Im_Lap_PAN = imfilter(I_PAN,fspecial('laplacian'));
Im_Lap_PAN = imfilter(I_PAN,fspecial('sobel'));
Im_Lap_PAN_rep = repmat(Im_Lap_PAN,[1,1,size(I_MS,3)]);

Im_Lap_F = zeros(size(I_F));
for idim=1:size(I_MS_LR,3)
    % Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('laplacian'));
    Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
end

% keyboard
%  Im_Lap_F(Im_Lap_F > 2^16) = 2^16;
%  Im_Lap_F(Im_Lap_F <- 2^16) = -2^16;

sCC=SCC(Im_Lap_PAN_rep,Im_Lap_F);

if nargout>=6
    Im_Lap_MS_LR=zeros(size(I_MS_LR));
    for idim=1:size(I_MS_LR,3)
        Im_Lap_MS_LR(:,:,idim)= imfilter(I_MS_LR(:,:,idim),fspecial('sobel'));
    end
    I_PAN_LR = LPF_filter(I_PAN,ratio,LR_method,true,GNyq_PAN);
    Im_Lap_PAN_LR = imfilter(I_PAN_LR,fspecial('sobel'));
    q=1;
    D_S_Sob=spatial_distortion_new(Im_Lap_F,Im_Lap_PAN,Im_Lap_MS_LR,Im_Lap_PAN_LR,q,ratio,Qblocks_size);
end
if nargout>=7
    Im_Lap_MS=zeros(size(I_MS));
    for idim=1:size(I_MS,3)
        Im_Lap_MS(:,:,idim)= imfilter(I_MS(:,:,idim),fspecial('sobel'));
    end
    I_PAN_F=LPF_filter(I_PAN,ratio,LR_method,false,GNyq_PAN);
    Im_Lap_PAN_F = imfilter(I_PAN_F,fspecial('sobel'));
    q=1;
    D_S_HR_Sob=spatial_distortion_new(Im_Lap_F,Im_Lap_PAN,Im_Lap_MS,Im_Lap_PAN_F,q,1,Qblocks_size);
end
if nargout>=8
    alpha=estimation_alpha(cat(3,ones(size(I_MS_LR,1),size(I_MS_LR,2)),I_MS_LR),I_PAN_LR,'global');
    I_Sim=sum(cat(3,ones(size(I_F,1),size(I_F,2)),I_F).*repmat(reshape(alpha,[1,1,length(alpha)]),[size(I_F,1),size(I_F,2),1]),3);
    Im_Lap_Sim=imfilter(I_Sim,fspecial('sobel'));
    sCC_Sim_Sob=SCC(Im_Lap_Sim,Im_Lap_PAN);
    % if nargout>=9
    %     sCC_Sim=SCC(I_Sim,I_PAN);
    % end
    if nargout>=9
        Im_Lap2_PAN = imfilter(I_PAN,fspecial('laplacian'));
        Im_Lap2_Sim = imfilter(I_Sim,fspecial('laplacian'));
        sCC_Sim_Lap = SCC(Im_Lap2_Sim,Im_Lap2_PAN);
    end
    if nargout>=10
        sCC_Sim_Sob_l=SCC(Im_Lap_Sim,Im_Lap_PAN,'local',Qblocks_size);
    end
    if nargout>=11
        UIQI_Sim_Sob_l=img_qi(Im_Lap_Sim,Im_Lap_PAN,Qblocks_size);
    end
end

% if nargout>=6
%     for idim=1:size(I_PAN,3)
%         Im_Lap_PAN(:,:,idim)= imfilter(I_PAN(:,:,idim),fspecial('sobel'));
%     end
%     for idim=1:size(I_MS_LR,3),
%         Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('sobel'));
%     end
%     sCC_Sobel=sum(sum(sum(Im_Lap_PAN.*Im_Lap_F)));
%     sCC_Sobel=sCC_Sobel/sqrt(sum(sum(sum((Im_Lap_PAN.^2)))));
%     sCC_Sobel=sCC_Sobel/sqrt(sum(sum(sum((Im_Lap_F.^2)))));
% end

end