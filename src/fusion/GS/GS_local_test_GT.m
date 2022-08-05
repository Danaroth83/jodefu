%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           GS fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Gram-Schmidt (GS) transformation.
% 
% Interface:
%           I_Fus_GS = GS(I_MS,I_PAN)
%
% Inputs:
%           I_MS:       MS image upsampled at PAN scale;
%           I_PAN:      PAN image.
%
% Outputs:
%           I_Fus_GS:   GS pasharpened image.
% 
% References:
%           [Laben00]   C. A. Laben and B. V. Brower, “Process for enhancing the spatial resolution of multispectral imagery using pan-sharpening,” Eastman
%                       Kodak Company, Tech. Rep. US Patent # 6,011,875, 2000.
%           [Aiazzi07]  B. Aiazzi, S. Baronti, and M. Selva, “Improving component substitution Pansharpening through multivariate regression of MS+Pan
%                       data,” IEEE Transactions on Geoscience and Remote Sensing, vol. 45, no. 10, pp. 3230–3239, October 2007.
%           [Vivone14]   G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                        IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_Fus_GS,g] = GS_local_test_GT(I_MS,I_PAN)

imageLR = double(I_MS);
imageHR = double(I_PAN);

%%% Remove means from imageLR
imageLR0 = zeros(size(I_MS));
for ii = 1 : size(I_MS,3), imageLR0(:,:,ii) = imageLR(:,:,ii) - mean2(imageLR(:,:,ii)); end

%%% Sintetic intensity
I = mean(imageLR,3); 

%% Coefficients
edge_extension=[floor((block_win-block_shift)/2);ceil((block_win-block_shift)/2)];
imageLR0_ext=edge_extend(imageLR0,edge_extension);
imageHR_ext=edge_extend(imageHR,edge_extension);
I_ext=edge_extend(I,edge_extension);
[L1e,L2e,Nb]=size(imageLR0_ext);
stepx=floor((L1e-block_win)/block_shift+1);
stepy=floor((L2e-block_win)/block_shift+1);
g = ones(stepx,stepy,size(I_MS,3)+1);
I_Fus_GSA=zeros(size(I_MS));
for ii = 1 : stepx
    iA=(ii-1)*block_shift+(1:block_win);
    iA1=(ii-1)*block_shift+(1:block_shift);
    iA2=ceil((block_win-block_shift)/2)+(1:length(iA1));
    for jj = 1 : stepy
        iB=(jj-1)*block_shift+(1:block_win);
        iB1=(jj-1)*block_shift+(1:block_shift);
        iB2=floor((block_win-block_shift)/2)+(1:length(iB1));
        Iblock = I_ext(iA,iB); 
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
        Vhat=V + repmat(delta,[1,1,size(I_MS,3)+1]) .* g(ii,jj,:);
        I_Fus_GSA(iA1,iB1,:)=Vhat(:,:,2:end);
    end
end

% keyboard
V_hat = V + deltam .* gm;

%%% Reshape fusion result
I_Fus_GS = reshape(V_hat(:,2:end),[size(I_MS,1) size(I_MS,2) size(I_MS,3)]);

% Final Mean Equalization
for ii = 1 : size(I_MS,3)
    h = I_Fus_GS(:,:,ii);
    I_Fus_GS(:,:,ii) = h - mean2(h) + mean2(imageLR(:,:,ii));
end

end