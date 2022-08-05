function y=idwt_custom(x,lo,hi,lev,sizes)

    if nargin<=3, lev=1; end
    if nargin<=4, sizes=size(x); end
    
    sizex=size(x,1);
    x=reshape(x,size(x),[]);
    
    Lfh=length(hi)-1;
    Lx=sizes(1);
    odd=zeros(1,lev);
    tail=zeros(1,lev+1);
    tail(lev+1)=size(x,1);

    
    for jj=lev:-1:1
        odd(jj)= (mod(Lx,2)==1);
        tail(jj)=tail(jj+1)-Lfh-floor(Lx/2);
        Lx=ceil(Lx/2);
    end
    
    y=[];
    for ii=1:size(x,2)
        ca=x(1:tail(jj),ii);
        for jj=1:lev
            cd=x(tail(jj)+1:tail(jj+1),ii);
            ca=idwt(ca,cd,lo,hi,'mode','per');
            if odd(jj), ca=ca(1:end-1); end
        end
        y=cat(2,y,ca);
    end
    
    sizex(1)=size(y,1);
    y=reshape(y,sizex);
    
end