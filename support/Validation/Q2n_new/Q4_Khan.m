function [ff,fx]=Q4_Khan(reference,fusa,sizeB,shift, extend_flag,norm_flag)
%
% usage: [ff,fx]=Q4(reference,fusa,sizeB,shift);
% ff: exit value of Q4 index (average value on sizexsize blocks). This is the main output of the function.
% fx: map of Q4 values computed on sizexsize blocks
% reference: reference 4-band image
% fusa: 4-band image to be evaluated
% shift: block displacement (in pixels)
% 	If shift=1, the fx map contains Q4 values for every image pixel
%	If shift=size, Q4 values are computed over adjacent, non-overlapped blocks
% extend_flag: calculate border values by periodization (different
%   behaviour between left and right sizes)
% normalization: normalize mean and variance w.r.t. reference signal (no
% simmetry between reference and fused)

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

if extend_flag
    
    stepx=ceil(N1/shift);
    stepy=ceil(N2/shift);
    
    
    if stepy<=0
        stepy=1;
        stepx=1;
    end;
    est1=(stepx-1)*shift+sizeB-N1;
    est2=(stepy-1)*shift+sizeB-N2;
    
    if sum([(est1~=0),(est2~=0)])>0
        N1_old=floor(est1/2);
        N2_old=floor(est2/2);
        ia1=zeros(N1+est1,N2+est2);
        ia1(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=ima1;
        ia1(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=ia1(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        ia1(:,N2_old:-1:1)=ia1(:,N2_old+1:N2_old+floor(est2/2));
        ia1(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=ia1(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        ia1(N1_old:-1:1,:)=ia1(N1_old+1:N1_old+floor(est1/2),:);
        ima1=ia1;
        
        ib1=zeros(N1+est1,N2+est2);
        ib1(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=imb1;
        ib1(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=ib1(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        ib1(:,N2_old:-1:1)=ib1(:,N2_old+1:N2_old+floor(est2/2));
        ib1(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=ib1(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        ib1(N1_old:-1:1,:)=ib1(N1_old+1:N1_old+floor(est1/2),:);
        imb1=ib1;
        
        
        ic1=zeros(N1+est1,N2+est2);
        ic1(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=imc1;
        ic1(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=ic1(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        ic1(:,N2_old:-1:1)=ic1(:,N2_old+1:N2_old+floor(est2/2));
        ic1(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=ic1(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        ic1(N1_old:-1:1,:)=ic1(N1_old+1:N1_old+floor(est1/2),:);
        imc1=ic1;
        
        id1=zeros(N1+est1,N2+est2);
        id1(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=imd1;
        id1(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=id1(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        id1(:,N2_old:-1:1)=id1(:,N2_old+1:N2_old+floor(est2/2));
        id1(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=id1(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        id1(N1_old:-1:1,:)=id1(N1_old+1:N1_old+floor(est1/2),:);
        imd1=id1;
        
        ia2=zeros(N1+est1,N2+est2);
        ia2(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=ima2;
        ia2(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=ia2(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        ia2(:,N2_old:-1:1)=ia2(:,N2_old+1:N2_old+floor(est2/2));
        ia2(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=ia2(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        ia2(N1_old:-1:1,:)=ia2(N1_old+1:N1_old+floor(est1/2),:);
        ima2=ia2;
        
        ib2=zeros(N1+est1,N2+est2);
        ib2(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=imb2;
        ib2(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=ib2(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        ib2(:,N2_old:-1:1)=ib2(:,N2_old+1:N2_old+floor(est2/2));
        ib2(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=ib2(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        ib2(N1_old:-1:1,:)=ib2(N1_old+1:N1_old+floor(est1/2),:);
        imb2=ib2;
        
        
        ic2=zeros(N1+est1,N2+est2);
        ic2(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=imc2;
        ic2(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=ic2(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        ic2(:,N2_old:-1:1)=ic2(:,N2_old+1:N2_old+floor(est2/2));
        ic2(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=ic2(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        ic2(N1_old:-1:1,:)=ic2(N1_old+1:N1_old+floor(est1/2),:);
        imc2=ic2;
        
        id2=zeros(N1+est1,N2+est2);
        id2(N1_old+1:N1_old+N1,N2_old+1:N2_old+N2)=imd2;
        id2(:,N2_old+N2+1:N2_old+N2+ceil(est2/2))=id2(:,N2_old+N2:-1:N2_old+N2-ceil(est2/2)+1);
        id2(:,N2_old:-1:1)=id2(:,N2_old+1:N2_old+floor(est2/2));
        id2(N1_old+N1+1:N1_old+N1+ceil(est1/2),:)=id2(N1_old+N1:-1:N1_old+N1-ceil(est1/2)+1,:);
        id2(N1_old:-1:1,:)=id2(N1_old+1:N1_old+floor(est1/2),:);
        imd2=id2;
    end
else
    %%% Without extension
    stepx=floor((N1-sizeB)/shift+1);
    stepy=floor((N2-sizeB)/shift+1);
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
        if norm_flag
            [o1,o2,o3,o4]=q_universal(rit_ima1,rit_imb1,rit_imc1,rit_imd1,rit_ima2,rit_imb2,rit_imc2,rit_imd2,sizeB,sizeB);
        else
            [o1,o2,o3,o4]=q_universal_nonorm(rit_ima1,rit_imb1,rit_imc1,rit_imd1,rit_ima2,rit_imb2,rit_imc2,rit_imd2,sizeB,sizeB);
        end
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