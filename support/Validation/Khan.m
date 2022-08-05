%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%          Khan indexes
%
% Interface:
%           [D_lambda_Khan,D_S_Khan]=Khan(I_F,I_MS,I_PAN,I_MS_LR,sensor,tag,ratio,Qblocks_size)
%
% Inputs:
%           I_F:                Pansharpened image;
%           I_MS:               MS image resampled to panchromatic scale;
%           I_MS_LR:            MS image;
%           I_PAN:              Panchromatic image;
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           tag:                Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value.
%           Q_blocks_size:      Block size of the Q-index locally applied;
%
% Outputs:
%           D_lambda_Khan:      D_lambda_Khan index;
%           D_s_Khan:           D_s_Khan index.
%
% References:
%           [Khan09]            Khan, L. Alparone, and J. Chanussot, “Pansharpening quality
%                               assessment using the modulation transfer functions of instruments,”
%                               IEEE Trans. Geosci. Remote Sens., vol. 11, no. 47, pp. 3880–3891,
%                               Nov. 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D_lambda_Khan,D_S_Khan]=Khan(I_F,I_MS,I_PAN,I_MS_LR,sensor,tag,ratio,Qblocks_size)


% keyboard
I_PAN_rep=repmat(I_PAN,[1 1 size(I_MS,3)]);
I_PAN_LR = imresize(ATWT_dec(I_PAN,ratio),1/ratio,'nearest');
I_PAN_LR_rep=repmat(I_PAN_LR,[1 1 size(I_MS,3)]);
% t = imresize(I_PAN_LR,1/ratio,'nearest');

%% MS_LP and MS_HP
% I_F_LP = MTF_noresize(I_F,sensor,tag,ratio);
I_F_LP = MTF(I_F,sensor,tag,ratio);
I_F_HP = I_F - I_F_LP;

%% MS_LR_LP and MS_LR_HP
% I_MS_LR_LP = MTF_noresize(I_MS_LR,sensor,tag,ratio);
I_MS_LR_LP = MTF(I_MS_LR,sensor,tag,ratio);
I_MS_LR_HP = I_MS_LR - I_MS_LR_LP;


%% PAN_LP and PAN_HP
%%%Khan version (almost ideal filters)
t = imresize(I_PAN,1/ratio,'nearest');
if (2^round(log2(ratio)) ~= ratio)
%     I_PAN_LP = interp23tapGeneral(t,ratio);
        I_PAN_LP = imresize(t,[size(I_PAN,1),size(I_PAN,2)]);
else
    I_PAN_LP = interp23tap(t,ratio);
end
I_PAN_LP_rep = repmat(I_PAN_LP,[1 1 size(I_MS,3)]);
I_PAN_HP_rep = I_PAN_rep - I_PAN_LP_rep;

t = imresize(I_PAN_LR,1/ratio,'nearest');
if (2^round(log2(ratio)) ~= ratio)
%     I_PAN_LR_LP = interp23tapGeneral(t,ratio);
        I_PAN_LR_LP = imresize(t,[size(I_PAN_LR,1),size(I_PAN_LR,2)]);
else
    I_PAN_LR_LP = interp23tap(t,ratio);
end
% keyboard
I_PAN_LR_LP_rep = repmat(I_PAN_LR_LP,[1 1 size(I_MS,3)]);
I_PAN_LR_HP_rep = I_PAN_LR_rep - I_PAN_LR_LP_rep;


% %%% Aiazzi version
% B. Aiazzi, L. Alparone, , S. Baronti, R. Carlà, A. Garzelli, and
% L. Santurri, “Full scale assessment of pansharpening methods and data
% products,” in SPIE Remote Sensing, 2014, pp. 924402–924402.

% % I_PAN_LP_rep  = MTF_noresize(I_PAN_rep,sensor,tag,ratio);
% I_PAN_LP_rep  = MTF(I_PAN_rep,sensor,tag,ratio);
% I_PAN_HP_rep = I_PAN_rep - I_PAN_LP_rep;
%
% I_PAN_LR_LP_rep  = MTF_noresize(I_PAN_LR_rep,sensor,tag,ratio);
% I_PAN_LR_HP_rep = I_PAN_LR_rep - I_PAN_LR_LP_rep;


cd ../Quality_Indices

%% D_lambda
I_F_LP_down = imresize(I_F_LP,1/ratio,'nearest');

D_lambda_Khan = 1 - q2n(I_F_LP_down,I_MS_LR,Qblocks_size,Qblocks_size);
%% Watch out: according to the original paper the first argument is the I_MS_LP (q2n is not symmetric)
% D_lambda_Khan2= 1 - q2n(I_MS_LR,I_F_LP,Qblocks_size,Qblocks_size);

%% D_S
Q_orig_HR = zeros(1,size(I_MS,3));
Q_orig_LR = zeros(1,size(I_MS,3));
D_S_Khan_Band=0;
for ii=1:size(I_MS,3),
    Q_orig_HR(ii) =img_qi(I_F_HP(:,:,ii),I_PAN_HP_rep(:,:,ii),Qblocks_size);
    Q_orig_LR(ii) = img_qi(I_MS_LR_HP(:,:,ii),I_PAN_LR_HP_rep(:,:,ii),ceil(Qblocks_size/ratio));
    D_S_Khan_Band(ii)=  abs(Q_orig_LR(ii)-Q_orig_HR(ii));  
end

D_S_Khan=mean(D_S_Khan_Band);

% D_S_Khan_HR = 1 - Q_HR;
% D_S_Khan_LR = 1 - Q_LR;
% D_S_Khan=(1-D_S_Khan_LR)-(1-D_S_Khan_HR);

cd ../Quality_Indices_New



