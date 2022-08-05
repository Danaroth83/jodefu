function psnr_index = PSNR ( I_GT , I_F, L )

if isempty(L)
    psnr_index=20*log10(max(cat(1,I_GT(:),I_F(:)))/sqrt(MSE(I_GT,I_F)));
else
    psnr_index=20*log10((2^L-1)/sqrt(MSE(I_GT,I_F)));
end

end