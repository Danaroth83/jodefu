function [ stdv, avg ] = std_block( I_in, B )
%MEANCOV_NORM_BLOCK Finds standard deviation of a block processed image

[L1,L2]=size(I_in);
edge=floor(B/2);
I_in=edge_extend(I_in,edge);
Bsq=B^2;
avg=zeros(L1,L2);
stdv=zeros(L1,L2);
for i1=1:L1
    iA=(i1-1)+(1:B);
    for i2=1:L2
        block=reshape(I_in(iA,(i2-1)+(1:B)),[Bsq,1]);
        avg(i1,i2)=sum(block);
        stdv(i1,i2)=sum(block.^2);
    end
end
avg=avg/Bsq;
stdv=sqrt((stdv-avg.^2)/(Bsq-1));

end

