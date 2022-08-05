function I_out=opdir_fusion(I_in,lpfilter,spectralweight)

I_PAN= permute(sum(I_in.*shiftdim(spectralweight,-2),3),[1,2,4,3]);

I_MS=zeros(size(I_in));
for ii=1:size(I_in,3)
    I_MS(:,:,ii) = imfilter(I_in(:,:,ii),lpfilter(:,:,ii),'replicate','conv');
end
I_out=cat(3,I_PAN,I_MS);

end