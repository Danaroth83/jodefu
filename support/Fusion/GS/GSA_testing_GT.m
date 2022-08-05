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
%           [Aiazzi07]   B. Aiazzi, S. Baronti, and M. Selva, “Improving component substitution Pansharpening through multivariate regression of MS+Pan
%                        data,” IEEE Transactions on Geoscience and Remote Sensing, vol. 45, no. 10, pp. 3230–3239, October 2007.
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_Fus_GSA,alpha,g] = GSA_testing_GT(I_GT,I_PAN,I_MS_LR,block_win,block_shift)

if nargin<=3, block_win=33; end
if nargin<=4, block_shift=1; end

imageLR = double(I_GT);
imageHR = double(I_PAN);
imageLR_LP = double(I_MS_LR);

%%% Remove means from imageLR
imageLR0 = zeros(size(I_GT));
for ii = 1 : size(I_GT,3), imageLR0(:,:,ii) = imageLR(:,:,ii) - mean2(imageLR(:,:,ii)); end

%%% Remove means from imageLR_LP
imageLR_LP0 = zeros(size(I_MS_LR));
for ii = 1 : size(I_MS_LR,3), imageLR_LP0(:,:,ii) = imageLR_LP(:,:,ii) - mean2(imageLR_LP(:,:,ii)); end


%%% Sintetic intensity
imageHR0 = imageHR - mean2(imageHR);
% imageHR0 = LowPass_Bspline(imageHR0,ratio);
% imageHR0 = LPfilterWithoutDec(imageHR0,ratio);
[~,alpha] = estimation_alpha(cat(3,imageLR,ones(size(I_GT,1),size(I_GT,2))),imageHR0,'local',block_win,block_shift);

%% Coefficients
edge_extension=[floor((block_win-block_shift)/2);ceil((block_win-block_shift)/2)];
imageLR0_ext=edge_extend(imageLR0,edge_extension);
imageHR_ext=edge_extend(imageHR,edge_extension);
[L1e,L2e,Nb]=size(imageLR0_ext);
stepx=floor((L1e-block_win)/block_shift+1);
stepy=floor((L2e-block_win)/block_shift+1);
g = ones(stepx,stepy,size(I_GT,3)+1);
I_Fus_GSA=zeros(size(I_GT));
iA2=ceil((block_win-block_shift)/2)+(1:block_shift);
iB2=floor((block_win-block_shift)/2)+(1:block_shift);
for ii = 1 : stepx
    iA=((ii-1)*block_shift)+(1:block_win);
    iA1=(ii-1)*block_shift+(1:block_shift);
    for jj = 1 : stepy
        iB=((jj-1)*block_shift)+(1:block_win);
        iB1=(jj-1)*block_shift+(1:block_shift);
        hblock=imageLR0_ext(iA,iB,:);
        Iblock = sum(cat(3,hblock,ones(block_win)) .* repmat(alpha(ii,jj,:),[block_win,block_win,1]),3); 
        I0block = Iblock - mean2(Iblock);
        imageLR0_min=imageLR0_ext(iA,iB,:);
        for zz = 1 : Nb
            h = imageLR0_min(:,:,zz);
            c = cov(I0block(:),h(:));
            g(ii,jj,zz+1) = c(1,2)/var(I0block(:));
        end
        imageHRblock=imageHR(iA1,iB1)-mean2(imageHR_ext(iA,iB));
        I0minblock=I0block(iA2,iB2);
        delta=imageHRblock-I0minblock;
        V=cat(3,I0minblock,imageLR0_min(iA2,iB2,:));
        Vhat=V + repmat(delta,[1,1,size(I_GT,3)+1]) .* g(ii,jj,:);
        I_Fus_GSA(iA1,iB1,:)=Vhat(:,:,2:end);
    end
end

% Final Mean Equalization
for ii = 1 : size(I_GT,3)
    h = I_Fus_GSA(:,:,ii);
    I_Fus_GSA(:,:,ii) = h - mean2(h) + mean2(imageLR(:,:,ii));
end

end