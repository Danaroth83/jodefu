%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%          Zhou index (QNR). 
% 
% Interface:
%           [Zhou_index,D_S_Zhou,D_lambda_Zhou] = Zhou(I_F,I_MS,I_PAN)
%
% Inputs:
%           I_F:                Pansharpened image;
%           I_MS:               MS image resampled to panchromatic scale;
%           I_PAN:              Panchromatic image;
% 
% Outputs:
%           Zhou_index:         Zhou index;
%           D_lambda_Zhou:      D_lambda_Zhou index;
%           D_s_Zhou:           D_s_Zhou index.
% 
% References:
%           [Zhou98]            J. Zhou, D. L. Civco, and J. A. Silander, “A wavelet transform method to merge landsat TM and SPOT panchromatic data,” Int. J. Remote Sens., vol. 19, no. 4, 
%                               pp. 743–757, May 1998
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D_S_Zhou,D_lambda_Zhou] = Zhou(I_F,I_MS,I_PAN)


%% DS

I_PAN = repmat(I_PAN,[1 1 size(I_MS,3)]);

Im_Lap_PAN = zeros(size(I_PAN));
for idim=1:size(I_PAN,3),
     Im_Lap_PAN(:,:,idim)= imfilter(I_PAN(:,:,idim),fspecial('laplacian'));
end

Im_Lap_F = zeros(size(I_F));
for idim=1:size(I_F,3),
    Im_Lap_F(:,:,idim)= imfilter(I_F(:,:,idim),fspecial('laplacian'));
end
% keyboard
sCC = SCC(Im_Lap_PAN,Im_Lap_F);

D_S_Zhou = 1-sCC;

%% D_lambda
err=abs(I_F-I_MS);

D_lambda_Zhou=mean(err(:));

