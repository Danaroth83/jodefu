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
function [D_lambda,D_S,QNR_index,Q_PAN_FR,Q_PAN_RR,SAM_index,sCC] = indexes_evaluation_FS_test(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,band_overlap_cell,LR_method,Qblocks_shift,Qblocks_shiftl,q)
if nargin<=10 || isempty(LR_method)
    LR_method='resize';
end
if nargin<=11, Qblocks_shift=Qblocks_size; end
if nargin<=12, Qblocks_shiftl=Qblocks_size; end
if nargin<=13, q=1; end

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

Nb_MS=size(I_PAN,3);
QNR_index=zeros(1,Nb_MS);
D_lambda=zeros(1,Nb_MS);
D_S=zeros(1,Nb_MS);
Q_PAN_FR=zeros(1,Nb_MS);
Q_PAN_RR=zeros(1,Nb_MS);
SAM_index=zeros(1,Nb_MS);
sCC=zeros(1,Nb_MS);
bo_len=zeros(1,Nb_MS);
for ii=1:Nb_MS
    I_PAN_temp=I_PAN(:,:,ii);
    I_MS_LR_temp=I_MS_LR(:,:,band_overlap_cell{ii});
    I_MS_temp=I_MS(:,:,band_overlap_cell{ii});
    I_F_temp=I_F(:,:,band_overlap_cell{ii});
    GNyq_temp=GNyq(band_overlap_cell{ii});
    GNyq_PAN_temp=GNyq_PAN(ii);
    [QNR_index(ii),D_lambda(ii),D_S(ii),~,Q_PAN_FR(ii),Q_PAN_RR(ii)]= QNR_mod(I_F_temp,I_MS_LR_temp,I_MS_temp,I_PAN_temp,GNyq_PAN_temp,ratio,Qblocks_size,LR_method,1,q,1,1,Qblocks_shift,Qblocks_shiftl);

    %     I_F_temp(I_F_temp > 2^L) = 2^L;
    %     I_F_temp(I_F_temp < 0) = 0;

    I_PAN_temp = repmat(I_PAN_temp,[1 1 size(I_MS_temp,3)]);

    % I_F_D = MTF(I_F,sensor,tag,ratio);
    I_F_D = MTF_filter(I_F_temp,GNyq_temp,ratio,41);
    I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

    SAM_index(ii) = SAM(I_MS_LR_temp,I_F_D); %toolbox

    sCC(ii)=SCC(I_PAN_temp,I_F_temp,'laplacian','global1');
    
    bo_len(ii)=length(band_overlap_cell{ii});
end
bands=length(cell2mat(band_overlap_cell));
QNR_index=sum(QNR_index.*bo_len)/bands;
D_lambda=sum(D_lambda.*bo_len)/bands;
D_S=sum(D_S.*bo_len)/bands;
Q_PAN_FR=sum(Q_PAN_FR.*bo_len)/bands;
Q_PAN_RR=sum(Q_PAN_RR.*bo_len)/bands;
SAM_index=sum(SAM_index.*bo_len)/bands;
sCC=sum(sCC.*bo_len)/bands;

end