function [ff,fx]=Q4(reference,fusa,sizeB,shift)
%
% usage: [ff,fx]=Q4(reference,fusa,sizeB,shift);
% ff: exit value of Q4 index (average value on sizexsize blocks). This is the main output of the function.
% fx: map of Q4 values computed on sizexsize blocks
% reference: reference 4-band image
% fusa: 4-band image to be evaluated
% shift: block displacement (in pixels)
% 	If shift=1, the fx map contains Q4 values for every image pixel  
%	If shift=size, Q4 values are computed over adjacent, non-overlapped blocks
%

N1 = size(reference,1);
N2 = size(reference,2);

reference=double(reference);
fusa=double(fusa);
ima1=squeeze(reference(:,:,1));
imb1=squeeze(reference(:,:,2));
imc1=squeeze(reference(:,:,3));
imd1=squeeze(reference(:,:,4));
ima2=squeeze(fusa(:,:,1));
imb2=squeeze(fusa(:,:,2));
imc2=squeeze(fusa(:,:,3));
imd2=squeeze(fusa(:,:,4));

%[N1,N2]=sizeB(ima1);
% stepx=round(((N1-sizeB)/shift)-1)+2;
% stepy=round(((N2-sizeB)/shift)-1)+2;

stepx=ceil(N1/shift);
stepy=ceil(N2/shift);

if stepy<=0
    stepy=1;
    stepx=1;
end;
valori=zeros(stepx,stepy);

% est1=stepx*sizeB-N1;
% est2=stepy*sizeB-N2;
est1=(stepx-1)*shift+sizeB-N1;
est2=(stepy-1)*shift+sizeB-N2;
if sum([(est1~=0),(est2~=0)])>0
  ia1=zeros(N1+est1,N2+est2);
  ia1(1:N1,1:N2)=ima1;
  ia1(:,N2+1:N2+est2)=ia1(:,N2:-1:N2-est2+1);
  ia1(N1+1:N1+est1,:)=ia1(N1:-1:N1-est1+1,:);
  ima1=ia1;
  
  ib1=zeros(N1+est1,N2+est2);
  ib1(1:N1,1:N2)=imb1;
  ib1(:,N2+1:N2+est2)=ib1(:,N2:-1:N2-est2+1);
  ib1(N1+1:N1+est1,:)=ib1(N1:-1:N1-est1+1,:);
  imb1=ib1;
  
  ic1=zeros(N1+est1,N2+est2);
  ic1(1:N1,1:N2)=imc1;
  ic1(:,N2+1:N2+est2)=ic1(:,N2:-1:N2-est2+1);
  ic1(N1+1:N1+est1,:)=ic1(N1:-1:N1-est1+1,:);
  imc1=ic1;
  
  id1=zeros(N1+est1,N2+est2);
  id1(1:N1,1:N2)=imd1;
  id1(:,N2+1:N2+est2)=id1(:,N2:-1:N2-est2+1);
  id1(N1+1:N1+est1,:)=id1(N1:-1:N1-est1+1,:);
  imd1=id1;
  
  ia2=zeros(N1+est1,N2+est2);
  ia2(1:N1,1:N2)=ima2;
  ia2(:,N2+1:N2+est2)=ia2(:,N2:-1:N2-est2+1);
  ia2(N1+1:N1+est1,:)=ia2(N1:-1:N1-est1+1,:);
  ima2=ia2;
  
  ib2=zeros(N1+est1,N2+est2);
  ib2(1:N1,1:N2)=imb2;
  ib2(:,N2+1:N2+est2)=ib2(:,N2:-1:N2-est2+1);
  ib2(N1+1:N1+est1,:)=ib2(N1:-1:N1-est1+1,:);
  imb2=ib2;
  
  ic2=zeros(N1+est1,N2+est2);
  ic2(1:N1,1:N2)=imc2;
  ic2(:,N2+1:N2+est2)=ic2(:,N2:-1:N2-est2+1);
  ic2(N1+1:N1+est1,:)=ic2(N1:-1:N1-est1+1,:);
  imc2=ic2;
  
  id2=zeros(N1+est1,N2+est2);
  id2(1:N1,1:N2)=imd2;
  id2(:,N2+1:N2+est2)=id2(:,N2:-1:N2-est2+1);
  id2(N1+1:N1+est1,:)=id2(N1:-1:N1-est1+1,:);
  imd2=id2; 
end


qs=zeros(stepx,stepy);
qi=zeros(stepx,stepy);
qj=zeros(stepx,stepy);
qk=zeros(stepx,stepy);

for j=1:stepx
    for i=1:stepy
        rit_ima1=ima1(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_imb1=imb1(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_imc1=imc1(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_imd1=imd1(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_ima2=ima2(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_imb2=imb2(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_imc2=imc2(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        rit_imd2=imd2(((j-1)*shift)+1:((j-1)*shift)+sizeB,((i-1)*shift)+1:((i-1)*shift)+sizeB);
        [o1,o2,o3,o4]=q_universal(rit_ima1,rit_imb1,rit_imc1,rit_imd1,rit_ima2,rit_imb2,rit_imc2,rit_imd2,sizeB,sizeB);
        qs(j,i)=o1;
        qi(j,i)=o2;
        qj(j,i)=o3;
        qk(j,i)=o4;
    end;
    %disp(j)
end; 

f=sqrt(qs.^2+qi.^2+qj.^2+qk.^2);

%figure;imshow(f,[]);colorbar

%%%%%

ff=mean2(f);
fx=f;