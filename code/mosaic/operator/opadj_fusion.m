function I_out=opadj_fusion(I_in,lpfilter,spectralweights)

Nb_PAN=size(spectralweights,2);
I_out1 = sum(permute(I_in(:,:,1:Nb_PAN),[1,2,4,3]).*shiftdim(spectralweights,-2),4);

I_out2=zeros(size(I_out1));
for ii=1:size(lpfilter,3)
    I_out2(:,:,ii) = imfilter(I_in(:,:,Nb_PAN+ii),lpfilter(:,:,ii),'corr','replicate');
end

I_out = I_out1+I_out2;

end