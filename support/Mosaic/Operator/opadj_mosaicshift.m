function x = opadj_mosaicshift(y,cfa,cfa_shift,padding)
    
    x=zeros(size(y,1),size(y,2),size(cfa,3));
    for ii=1:size(cfa,3)
        x(:,:,ii)=circshift(y,-cfa_shift(ii,:));
    end
    x=x(padding(1)+1:end-padding(3),padding(2)+1:end-padding(4),:);
    x=x.*cfa;
    
end