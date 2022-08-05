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
function [D_lambda,D_S,QNR_index,SAM_index,sCC,D_S_Sob,D_S_HR_Sob,sCC_Sim_Sob,sCC_Sim_Lap,sCC_Sim_Sob_l,UIQI_Sim_Sob_l] = indexes_evaluation_HS_test_grp(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,band_overlap_cell,LR_method)
if nargin<=11 || isempty(LR_method)
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

D_S=0; D_lambda=0; QNR_index=0;
for ii=1:length(band_overlap_cell)
    % [QNR_index_temp,D_lambda_temp,D_S_temp]= QNR(I_F,I_MS,I_PAN,sensor,ratio,Qblocks_size);
    [QNR_index_temp,D_lambda_temp,D_S_temp]= QNR_mod(I_F(:,:,band_overlap_cell{ii}),I_MS_LR(:,:,band_overlap_cell{ii}),I_MS(:,:,band_overlap_cell{ii}),I_PAN(:,:,ii),GNyq_PAN(ii),ratio,Qblocks_size,LR_method);
    QNR_index=QNR_index+QNR_index_temp*length(band_overlap_cell{ii});
    D_lambda=D_lambda+D_lambda_temp*length(band_overlap_cell{ii});
    D_S=D_S+D_S_temp*length(band_overlap_cell{ii});
end
QNR_index=QNR_index/length(cell2mat(band_overlap_cell));
D_lambda=D_lambda/length(cell2mat(band_overlap_cell));
D_S=D_S/length(cell2mat(band_overlap_cell));

%     I_F(I_F > 2^L) = 2^L;
%     I_F(I_F < 0) = 0;


% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

SAM_index=0;
for ii=1:length(band_overlap_cell)
    SAM_index_temp = SAM(I_MS_LR(:,:,band_overlap_cell{ii}),I_F_D(:,:,band_overlap_cell{ii})); %toolbox
    SAM_index=SAM_index+SAM_index_temp*length(band_overlap_cell{ii});
end
SAM_index=SAM_index/length(cell2mat(band_overlap_cell));

Im_Lap_PAN=zeros(size(I_PAN));
Im_Lap_PAN_rep=zeros(size(I_MS));
Im_Lap2_PAN=zeros(size(I_PAN));
for ii=1:size(I_PAN,3)
    Im_Lap_PAN(:,:,ii) = imfilter(I_PAN(:,:,ii),fspecial('sobel'));
    Im_Lap_PAN_rep(:,:,band_overlap_cell{ii}) = repmat(Im_Lap_PAN(:,:,ii),[1,1,length(band_overlap_cell{ii})]);
    if nargout>=10
        Im_Lap2_PAN(:,:,ii) = imfilter(I_PAN(:,:,ii),fspecial('laplacian'));
    end
end

Im_Lap_F = zeros(size(I_F));
for ii=1:size(I_MS_LR,3)
    Im_Lap_F(:,:,ii)= imfilter(I_F(:,:,ii),fspecial('sobel'));
end

% keyboard
%  Im_Lap_F(Im_Lap_F > 2^16) = 2^16;
%  Im_Lap_F(Im_Lap_F <- 2^16) = -2^16;
sCC=0;
for ii=1:length(band_overlap_cell)
    sCC_temp=SCC(Im_Lap_PAN_rep(:,:,band_overlap_cell{ii}),Im_Lap_F(:,:,band_overlap_cell{ii}));
    sCC=sCC+sCC_temp*length(band_overlap_cell{ii});
end
sCC=sCC/length(cell2mat(band_overlap_cell));

if nargout>=6
    Im_Lap_MS_LR=zeros(size(I_MS_LR));
    for ii=1:size(I_MS_LR,3)
        Im_Lap_MS_LR(:,:,ii)= imfilter(I_MS_LR(:,:,ii),fspecial('sobel'));
    end
    I_PAN_LR=zeros(size(I_PAN,1)/ratio,size(I_PAN,2)/ratio,size(I_PAN,3));
    Im_Lap_PAN_LR=zeros(size(I_PAN_LR));
    for ii=1:size(I_PAN,3)
        I_PAN_LR(:,:,ii) = LPF_filter(I_PAN(:,:,ii),ratio,LR_method,true,GNyq_PAN(ii));
        Im_Lap_PAN_LR(:,:,ii) = imfilter(I_PAN_LR(:,:,ii),fspecial('sobel'));            
    end
    q=1;
    D_S_Sob=0;
    for ii=1:length(band_overlap_cell)
        D_S_Sob_temp=spatial_distortion_new(Im_Lap_F(:,:,band_overlap_cell{ii}),Im_Lap_PAN(:,:,ii),Im_Lap_MS_LR(:,:,band_overlap_cell{ii}),Im_Lap_PAN_LR(:,:,ii),q,ratio,Qblocks_size);
        D_S_Sob=D_S_Sob+D_S_Sob_temp*length(band_overlap_cell{ii});
    end
    D_S_Sob=D_S_Sob/length(cell2mat(band_overlap_cell));
end
if nargout>=7
    Im_Lap_MS=zeros(size(I_MS));
    for ii=1:size(I_MS,3)
        Im_Lap_MS(:,:,ii)= imfilter(I_MS(:,:,ii),fspecial('sobel'));
    end
    I_PAN_F=zeros(size(I_PAN));
    Im_Lap_PAN_F=zeros(size(I_PAN));
    for ii=1:size(I_PAN,3)
        I_PAN_F(:,:,ii)=LPF_filter(I_PAN(:,:,ii),ratio,LR_method,false,GNyq_PAN(ii));
        Im_Lap_PAN_F(:,:,ii) = imfilter(I_PAN_F(:,:,ii),fspecial('sobel'));
    end
    q=1;
    D_S_HR_Sob=0;
    for ii=1:length(band_overlap_cell)
        D_S_HR_Sob_temp=spatial_distortion_new(Im_Lap_F(:,:,band_overlap_cell{ii}),Im_Lap_PAN(:,:,ii),Im_Lap_MS(:,:,band_overlap_cell{ii}),Im_Lap_PAN_F(:,:,ii),q,1,Qblocks_size);
        D_S_HR_Sob=D_S_HR_Sob+D_S_HR_Sob_temp*length(band_overlap_cell{ii});
    end
    D_S_HR_Sob=D_S_HR_Sob/length(cell2mat(band_overlap_cell));
end
if nargout>=8
    sCC_Sim_Sob=0;
    % sCC_Sim=0;
    sCC_Sim_Lap=0;
    sCC_Sim_Sob_l=0;
    UIQI_Sim_Sob_l=0;
    for ii=1:length(band_overlap_cell)
        alpha=estimation_alpha(cat(3,ones(size(I_MS_LR,1),size(I_MS_LR,2)),I_MS_LR(:,:,band_overlap_cell{ii})),I_PAN_LR(:,:,ii),'global');
        I_Sim=sum(cat(3,ones(size(I_F,1),size(I_F,2)),I_F(:,:,band_overlap_cell{ii})).*repmat(reshape(alpha,[1,1,length(alpha)]),[size(I_F,1),size(I_F,2),1]),3);
        Im_Lap_Sim=imfilter(I_Sim,fspecial('sobel'));
        sCC_Sim_Sob_temp=SCC(Im_Lap_Sim,Im_Lap_PAN(:,:,ii));
        sCC_Sim_Sob=sCC_Sim_Sob+sCC_Sim_Sob_temp*length(band_overlap_cell{ii});
        % if nargout>=9
        %     sCC_Sim_temp=SCC(I_Sim,I_PAN(:,:,ii));
        %     sCC_Sim=sCC_Sim+sCC_Sim_temp*length(band_overlap_cell{ii});
        % end
        if nargout>=9
            Im_Lap2_Sim=imfilter(I_Sim,fspecial('laplacian'));
            sCC_Sim_Lap_temp=SCC(Im_Lap2_Sim,Im_Lap2_PAN(:,:,ii));
            sCC_Sim_Lap=sCC_Sim_Lap+sCC_Sim_Lap_temp*length(band_overlap_cell{ii});
        end
        if nargout>=10
            sCC_Sim_Sob_l_temp=SCC(Im_Lap_Sim,Im_Lap_PAN(:,:,ii),'local',Qblocks_size);
            sCC_Sim_Sob_l=sCC_Sim_Sob_l+sCC_Sim_Sob_l_temp*length(band_overlap_cell{ii});
        end
        if nargout>=11
            UIQI_Sim_Sob_l_temp=img_qi(Im_Lap_Sim,Im_Lap_PAN(:,:,ii),Qblocks_size);
            UIQI_Sim_Sob_l=UIQI_Sim_Sob_l+UIQI_Sim_Sob_l_temp*length(band_overlap_cell{ii});
        end
    end
    sCC_Sim_Sob=sCC_Sim_Sob/length(cell2mat(band_overlap_cell));
    % sCC_Sim=sCC_Sim/length(cell2mat(band_overlap_cell));
    sCC_Sim_Lap=sCC_Sim_Lap/length(cell2mat(band_overlap_cell));
    sCC_Sim_Sob_l=sCC_Sim_Sob_l/length(cell2mat(band_overlap_cell));
    UIQI_Sim_Sob_l=UIQI_Sim_Sob_l/length(cell2mat(band_overlap_cell));
end


end