function I_out=opadj_fusionfft(I_in,conj_fft2_KerBlu,spectralweight)

Nb_HR=size(spectralweight,2);
I_out1 = sum(permute(I_in(:,:,1:Nb_HR),[1,2,4,3]).*shiftdim(spectralweight,-2),4);
I_out2 = real(ifft2(fft2(I_in(:,:,Nb_HR+1:end)).*conj_fft2_KerBlu));
I_out = I_out1+I_out2;

end