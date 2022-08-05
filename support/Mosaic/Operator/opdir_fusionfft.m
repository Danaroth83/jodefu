function I_out=opdir_fusionfft(I_in,fft2_KerBlu,spectralweights)

I_PAN= permute(sum(I_in.*shiftdim(spectralweights,-2),3),[1,2,4,3]);
I_MS = real(ifft2(fft2_KerBlu.*fft2(I_in)));
I_out=cat(3,I_PAN,I_MS);

end