%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Quality with No Reference (QNR). 
% 
% Interface:
%           [QNR_index,D_lambda_index,D_s_index] = QNR(I_F,I_MS,I_PAN,sensor,ratio,S,p,q,alpha,beta)
%
% Inputs:
%           I_F:                Pansharpened image;
%           I_MS:               MS image resampled to panchromatic scale;
%           I_PAN:              Panchromatic image;
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value;
%           S:                  Block size (optional); Default value: 32;
%           p, q, alpha, beta:  Exponent values (optional); Default values: p = q = alpha = beta = 1.
% 
% Outputs:
%           QNR_index:          QNR index;
%           D_lambda_index:     D_lambda index;
%           D_s_index:          D_s index.
% 
% References:
%           [Alparone08]        L. Alparone, B. Aiazzi, S. Baronti, A. Garzelli, F. Nencini, and M. Selva, "Multispectral and panchromatic data fusion assessment without reference,"
%                               Photogrammetric Engineering and Remote Sensing, vol. 74, no. 2, pp. 193–200, February 2008. 
%           [Vivone14]          G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                               IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QNR_index,D_lambda_index,D_s_index,D_s_index2,Q_PAN_FR,Q_PAN_RR] = QNR_mod(I_F,I_MS_LR,I_MS,I_PAN,GNyq_PAN,ratio,S,LR_method,p,q,alpha,beta,shift,lshift)

if nargin <= 11, beta=1; end
if nargin <= 10, alpha=1; end
if nargin <= 9, q=1; end
if nargin <= 8, p=1; end
if nargin <= 7 || isempty(LR_method), LR_method='resize'; end
if nargin <= 6, S=32; end
if nargin <=12, shift=S; end
if nargin <=13, lshift=S; end

D_lambda_index = D_lambda_new(I_F,I_MS,S,p,lshift);

if nargout<=4
    % D_s_index = D_s_mod(I_F,I_MS,I_PAN,GNyq_PAN,ratio,S,q,LR_method); % toolbox
    D_s_index=DS_new(I_MS,I_F,I_PAN,ratio,GNyq_PAN,S,LR_method,q,shift);
else
    [D_s_index,Q_PAN_FR,Q_PAN_RR] = DS_new(I_MS,I_F,I_PAN,ratio,GNyq_PAN,S,LR_method,q,shift); % toolbox
end
if nargout==4
    D_s_index2 = DS_mod(I_MS_LR,I_F,I_PAN,ratio,GNyq_PAN,S,LR_method,q,shift); % new
else
    D_s_index2 = NaN;
end

QNR_index = (1-D_lambda_index)^alpha * (1-D_s_index)^beta;

end