%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           GSA fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Gram-Schmidt Adaptive (GSA) algorithm.
% 
% Interface:
%           I_Fus_GSA = GSA(I_MS,I_PAN,I_MS_LR,ratio)
%
% Inputs:
%           I_MS:       MS image upsampled at PAN scale;
%           I_PAN:      PAN image;
%           I_MS_LR:    MS image;
%           ratio:      Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_GSA:  GSA pasharpened image.
% 
% References:
%           [Aiazzi07]   B. Aiazzi, S. Baronti, and M. Selva, �Improving component substitution Pansharpening through multivariate regression of MS+Pan
%                        data,� IEEE Transactions on Geoscience and Remote Sensing, vol. 45, no. 10, pp. 3230�3239, October 2007.
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, �A Critical Comparison Among Pansharpening Algorithms�, 
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_Fus_GSA,alpha,g,Correlation_map] = GSA_local_CBD_framing(I_MS,I_PAN,~,ratio,block_win,block_shift,flag_downsample,I_PAN2)

if nargin<=4 || isempty(block_win), block_win=33; end
if nargin<=5 || isempty(block_shift), block_shift=1; end
if nargin<=6, flag_downsample=0; end

imageLR = double(I_MS);
imageHR = double(I_PAN);
imageHR_2 = double(I_PAN2);

%%% Remove means from imageLR
imageLR0 = zeros(size(I_MS));
for ii = 1 : size(I_MS,3), imageLR0(:,:,ii) = imageLR(:,:,ii) - mean2(imageLR(:,:,ii)); end

%%% Sintetic intensity
imageHR0 = imageHR - mean2(imageHR);
imageHR0_2 = imageHR_2 - mean2(imageHR_2);
if flag_downsample==1
    imageHR0 = LowPass_Bspline(imageHR0,ratio);
    imageHR0_2 = LowPass_Bspline(imageHR0_2,ratio);
else
    imageHR0 = LPfilterWithoutDec(imageHR0,ratio);
    imageHR0_2 = LPfilterWithoutDec(imageHR0_2,ratio);
end
[~,alpha] = estimation_alpha(cat(3,imageLR0,ones(size(I_MS,1),size(I_MS,2))),imageHR0,'local',block_win,block_shift);
[~,alpha_2] = estimation_alpha(cat(3,imageLR0,ones(size(I_MS,1),size(I_MS,2))),imageHR0_2,'local',block_win,block_shift);

%% Coefficients
edge_extension=[floor((block_win-block_shift)/2);ceil((block_win-block_shift)/2)];
imageLR0_ext=edge_extend(imageLR0,edge_extension);
imageHR_ext=edge_extend(imageHR,edge_extension);
imageHR0_ext=edge_extend(imageHR0,edge_extension);
imageHR_2_ext=edge_extend(imageHR_2,edge_extension);
imageHR0_2_ext=edge_extend(imageHR0_2,edge_extension);
[L1e,L2e,Nb]=size(imageLR0_ext);
stepx=floor((L1e-block_win)/block_shift+1);
stepy=floor((L2e-block_win)/block_shift+1);
g = ones(stepx,stepy,size(I_MS,3)+1);
I_Fus_GSA=zeros(size(I_MS));
Correlation_map=zeros(stepx,stepy);
Correlation_map_2=zeros(stepx,stepy);
for ii = 1 : stepx
    iA=(ii-1)*block_shift+(1:block_win);
    if ii==stepx
        iA1=((stepx-1)*block_shift+1):size(alpha,1);
    else
        iA1=(ii-1)*block_shift+(1:block_shift);
    end
    iA2=ceil((block_win-block_shift)/2)+(1:length(iA1));
    for jj = 1 : stepy
        iB=(jj-1)*block_shift+(1:block_win);
        if jj==stepy
            iB1=((stepy-1)*block_shift+1):size(alpha,2);
        else
            iB1=(jj-1)*block_shift+(1:block_shift);
        end
        iB2=floor((block_win-block_shift)/2)+(1:length(iB1));
        hblock=imageLR0_ext(iA,iB,:);
        Iblock = sum(cat(3,hblock,ones(block_win)) .* repmat(alpha(ii,jj,:),[block_win,block_win,1]),3); 
        I0block = Iblock - mean2(Iblock);
        I0minblock=I0block(iA2,iB2);
        Iblock_2 = sum(cat(3,hblock,ones(block_win)) .* repmat(alpha_2(ii,jj,:),[block_win,block_win,1]),3); 
        I0block_2 = Iblock_2 - mean2(Iblock_2);
        I0minblock_2=I0block_2(iA2,iB2);
        imageLR0_min=imageLR0_ext(iA,iB,:);
        imageHRblock=imageHR(iA1,iB1)-mean2(imageHR_ext(iA,iB));
        imageHRblock_2=imageHR_2(iA1,iB1)-mean2(imageHR_2_ext(iA,iB));
        Correlation_map(ii,jj)=corr2(I0block,imageHR0_ext(iA,iB));
        Correlation_map_2(ii,jj)=corr2(I0block_2,imageHR0_2_ext(iA,iB));
        if Correlation_map(ii,jj)>1.1*Correlation_map_2(ii,jj)
            for zz = 1 : Nb
                h = imageLR0_min(:,:,zz);
                c = cov(I0block(:),h(:));
                g(ii,jj,zz+1) = c(1,2)/var(I0block(:));
            end
            delta=imageHRblock-I0minblock;
            V=cat(3,I0minblock,imageLR0_min(iA2,iB2,:));
        else
            for zz = 1 : Nb
                h = imageLR0_min(:,:,zz);
                c = cov(I0block_2(:),h(:));
                g(ii,jj,zz+1) = c(1,2)/var(I0block_2(:));
            end
            delta=imageHRblock_2-I0minblock_2;
            V=cat(3,I0minblock_2,imageLR0_min(iA2,iB2,:));
        end
        Vhat=V + repmat(delta,[1,1,size(I_MS,3)+1]) .* g(ii,jj,:);
        I_Fus_GSA(iA1,iB1,:)=Vhat(:,:,2:end);
    end
end

% Final Mean Equalization
for ii = 1 : size(I_MS,3)
    h = I_Fus_GSA(:,:,ii);
    I_Fus_GSA(:,:,ii) = h - mean2(h) + mean2(imageLR(:,:,ii));
end

end