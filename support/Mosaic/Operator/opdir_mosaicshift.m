function x = opdir_mosaicshift(x,cfa,cfa_shift,padding)
    
    x=x.*cfa;
    x=padarray(x,padding(1:2),0,'pre');
    x=padarray(x,padding(3:4),0,'post');
    for ii=1:size(x,3)
        x(:,:,ii)=circshift(x(:,:,ii),cfa_shift(ii,:));
    end
    x=sum(x,3);
    
end