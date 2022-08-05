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
function [D_lambda,D_S,QNR_index,Q_PAN_FR,Q_PAN_RR,SAM_index,sCC] = indexes_evaluation_FS_test(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,GNyq,GNyq_PAN,ratio,Qblocks_size,LR_method,Qblocks_shift,Qblocks_shiftl,q)
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

[QNR_index,D_lambda,D_S,~,Q_PAN_FR,Q_PAN_RR]= QNR_mod(I_F,I_MS_LR,I_MS,I_PAN,GNyq_PAN,ratio,Qblocks_size,LR_method,1,q,1,1,Qblocks_shift,Qblocks_shiftl);

%     I_F(I_F > 2^L) = 2^L;
%     I_F(I_F < 0) = 0;

I_PAN = repmat(I_PAN,[1 1 size(I_MS,3)]);

% I_F_D = MTF(I_F,sensor,tag,ratio);
I_F_D = MTF_filter(I_F,GNyq,ratio,41);
I_F_D = I_F_D(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);

SAM_index = SAM(I_MS_LR,I_F_D); %toolbox

sCC=SCC(I_PAN,I_F,'laplacian','global1');

end