function Q = img_qi_global( I_F,I_PAN )
%IMG_QI_GLOBAL computes global UIQI index (I_PAN must be a 2D matrix)

[L1,L2,Nb]=size(I_F);
Dim=L1*L2;
I_1res=reshape(I_F,[Dim,Nb]);
I_2res=I_PAN(:);
m1=sum(I_1res,1);
m2=repmat(sum(I_2res,1),[1,Nb]);
s1=sum(I_1res.^2,1);
s2=repmat(sum(I_2res.^2,1),[1,Nb]);
s12=sum(I_1res.*repmat(I_2res,[1,Nb]),1);
mp=m1.*m2;
msq=m1.^2+m2.^2;
Q = 4*(Dim*s12-mp).*mp./(Dim*(s1+s2)-msq)./msq;

end

