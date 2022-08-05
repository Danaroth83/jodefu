%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Estimation coefficients linear regression model. 
% 
% Interface:
%           I_Fus_BDSD = BDSD(I_MS,I_PAN,ratio,S,sensor)
%
% Inputs:
%           I_MS:               MS image upsampled at PAN scale;
%           I_PAN:              PAN image;
%           type_estimation:    Type of estimation (i.e. local or global).
%
% Outputs:
%           alpha:              Coefficients estimated by the linear regression model.
% 
% References:
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha,alphas] = estimation_alpha(I_MS,I_PAN,type_estimation,block_win,block_shift)
if nargin<=2, type_estimation='global'; end
if nargin<=3, block_win=32; end
if nargin<=4, block_shift=block_win; end

if strcmp(type_estimation,'global')
    %%%%%%%% Global estimation
    IHc = reshape(I_PAN,[numel(I_PAN) 1]);
    ILRc = reshape(I_MS,[size(I_MS,1)*size(I_MS,2) size(I_MS,3)]);
    alpha = ILRc\IHc;
    alphas=alpha;
else
    %%%%%%%% Local estimation
    
    I_PAN=edge_extend(I_PAN,[ceil((block_win-block_shift)/2);floor((block_win-block_shift)/2)],'n');
    I_MS=edge_extend(I_MS,[ceil((block_win-block_shift)/2);floor((block_win-block_shift)/2)],'n');
    [L1e,L2e,Nb]=size(I_MS);
    stepx=floor((L1e-block_win)/block_shift+1);
    stepy=floor((L2e-block_win)/block_shift+1);
    block_win_sq=block_win^2;
    alphas = zeros(stepx,stepy,Nb);
    for ii = 1 : stepx
        iA=((ii-1)*block_shift)+(1:block_win);
        for jj = 1 : stepy
                iB=((jj-1)*block_shift)+(1:block_win);
                imHRbl = I_PAN(iA,iB);
                imageLRbl = I_MS(iA,iB,:);
                imageHRc = reshape(imHRbl,[block_win_sq 1]);
                ILRc = reshape(imageLRbl,[block_win_sq,Nb]);
                alphas(ii,jj,:) = ILRc\imageHRc;
        end
    end
    alpha = mean(reshape(alphas,[stepx*stepy,Nb]),1);        
end

end