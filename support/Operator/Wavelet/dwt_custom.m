function y=dwt_custom(x,lo,hi,lev)

    sizex=size(x);
    x=reshape(x,size(x,1),[]);

    y=[];
    for ii=1:size(x,2)
        vector=[];
        ca=x(:,ii);
        for jj=1:lev
            [ca,cd]=dwt(ca,lo,hi);
            vector=cat(1,cd,vector);
        end
        vector=cat(1,ca,vector);
        y=cat(2,y,vector);
    end
        
    sizex(1)=size(y,1);
    y=reshape(y,sizex);
    
end