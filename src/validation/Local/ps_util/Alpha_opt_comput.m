function [Alpha_opt,Alpha_opt_eq,Alpha_opt_eq2]=Alpha_opt_comput(I_MS,I_PAN,I_GT,ratio)


I_LP = genPAN_LP(I_PAN, 'ATWT', ratio);
Alpha_opt=(I_GT-I_MS)./(repmat((I_PAN - I_LP),[1 1 size(I_MS,3)]));

[Height,Width,Bands] = size(I_MS);

imageHR = (I_PAN - mean2(I_PAN))*...
        std2(I_LP)/std2(I_PAN) + mean2(I_LP);

Alpha_opt_eq=(I_GT-I_MS)./(repmat((imageHR - I_LP),[1 1 size(I_MS,3)]));

% %%% Equalization
imageHR = repmat(I_PAN,[1 1 size(I_MS,3)]);

for ii = 1 : size(I_MS,3)  
 imageHR(:,:,ii) = (imageHR(:,:,ii) - mean2(imageHR(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(imageHR(:,:,ii))) + mean2(I_MS(:,:,ii));  
end

Alpha_opt_eq2=(I_GT-I_MS)./(imageHR - repmat(I_LP,[1 1 size(I_MS,3)]));

