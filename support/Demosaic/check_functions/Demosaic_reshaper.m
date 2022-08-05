function I_out=Demosaic_reshaper(mask,I_in,options)
    
    perm=options.perm;
    
    I_out=zeros([size(I_in,1),size(I_in,2),max(perm(:))]);
    for kk=1:max(perm(:))
        I_out(:,:,kk)=mean(I_in(:,:,perm==kk),3);
    end
    for kk=1:size(mask,3)
        mask_idxnonzero=mask(:,:,kk)>0;
        I_out_temp=I_out(:,:,perm(kk));
        I_in_temp=I_in(:,:,kk);
        I_out_temp(mask_idxnonzero)=I_in_temp(mask_idxnonzero);
        I_out(:,:,perm(kk))=I_out_temp;
    end
    
    flip=options.flip;
    if flip(1)==1, I_out=I_out(end:-1:1,:,:); end
    if flip(2)==1, I_out=I_out(:,end:-1:1,:); end
    
    edge=options.edge;
    if edge(1)==1, I_out=I_out(1:end-1,:,:); end
    if edge(2)==1, I_out=I_out(:,1:end-1,:); end

end