function [ sc ] = hyp2n3D_scalar( v1, v2 )
%HYP2N3D_SCALAR Calculates the scalar product between 2^n-ons along the
%third dimension
%   v1 and v2 are 3d matrices whose third dimension is a power of 2

L=size(v1,3);
hL=L/2;
if L>2
    a1=v1(:,:,1);     a2=v1(:,:,2:hL);
    b1=v1(:,:,hL+1);  b2=-v1(:,:,hL+2:L);
    c1=v2(:,:,1);     c2=-v2(:,:,2:hL);
    d1=-v2(:,:,hL+1); d2=v2(:,:,hL+2:L);
    sc(:,:,1:hL)=hyp2n3D_scalar(cat(3,a1,a2),cat(3,c1,-c2))-...
        hyp2n3D_scalar(cat(3,d1,d2),cat(3,b1,b2));
    temp=hyp2n3D_scalar(cat(3,a1,-a2),cat(3,d1,-d2))+...
        hyp2n3D_scalar(cat(3,c1,c2),cat(3,b1,-b2));
    sc(:,:,hL+1:L)=cat(3,temp(:,:,1),-temp(:,:,2:hL));
elseif L==2
    sc=cat(3,sum(v1.*v2,3),sum(v1.*cat(3,-v2(:,:,2),v2(:,:,1)),3));
else
    sc=v1.*v2;
end
    
end
