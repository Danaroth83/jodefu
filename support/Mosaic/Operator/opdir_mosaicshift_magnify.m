function x = opdir_mosaicshift_magnify(x,cfa,cfa_shift,padding,filter)
    
    for ii=1:size(x,3)
        x(:,:,ii)=imfilter(x(:,:,ii),filter(:,:,ii),'circular','conv');
    end
    x=x.*cfa;
    x=padarray(x,padding(1:2),0,'pre');
    x=padarray(x,padding(3:4),0,'post');
    for ii=1:size(x,3)
        x(:,:,ii)=circshift(x(:,:,ii),cfa_shift(ii,:));
    end
    x=sum(x,3);
    
end