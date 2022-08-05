function [ Q2n_index, Q2n_map ] = Q2n_norm( I_1, I_2, B, K )
%Q2n Calculates the Q2^n quality index
%   I_1 is the pansharpened image set
%   I_2 are same-sized sets of reference images
%   B is the size of the each block used for the block statistic
%       (default B=32)
%   K is the amount of non overlapping pixels for each block shift
%       (default K=1)

%Default inputs
if nargin<4
    K=1;
end
if nargin<3
    B=32;
end

%[L1,L2,Nb]=size(I_1);

%Edge flipping
%I_1=edge_extend(I_1,[floor((B-K)/2);ceil((B-K)/2)],'n');
%I_2=edge_extend(I_2,[floor((B-K)/2);ceil((B-K)/2)],'n');
I_1=edge_extend(I_1,[floor(B-K);0],'r');
I_2=edge_extend(I_2,[floor(B-K);0],'r');
[L1e,L2e,Nb]=size(I_1);

%Zero-padding
Nbu=2^ceil(log2(Nb));
I_1=double(cat(3,I_1,zeros(L1e,L2e,Nbu-Nb)));
I_2=double(cat(3,I_2,zeros(L1e,L2e,Nbu-Nb)));

%Inizialization of block variables
L1f=floor((L1e-B)/K+1);
L2f=floor((L2e-B)/K+1);
m_z1=zeros(L1f,L2f,Nbu);
m_z2=zeros(L1f,L2f,Nbu);
corr_z1z2=zeros(L1f,L2f,Nbu);
corr_z1=zeros(L1f,L2f);
corr_z2=zeros(L1f,L2f);
Bsq=B^2;

%Block processing
fprintf('Block Processing for Q2^n index: %3d%%', 0);
for i1=1:L1f
    iA=(i1-1)*K+(1:B);
    for i2=1:L2f
        iB=(i2-1)*K+(1:B);
        
        I_1block=I_1(iA,iB,:);
        I_2block=I_2(iA,iB,:);
        
        %Normalization
        m1_s=sum(sum(I_1block,1),2)/Bsq;
        m1_s_mat=repmat(reshape(m1_s,[1,1,Nbu]),[B,B,1]);
        cov1_s=sum(sum((I_1block-m1_s_mat).^2,1),2)/Bsq;
        cov1_s(cov1_s==0)=1;
        cov1_s_mat=repmat(reshape(cov1_s,[1,1,Nbu]),[B,B,1]);
        I_1block=(I_1block-m1_s_mat)./cov1_s_mat+1;
        I_2block=(I_2block-m1_s_mat)./cov1_s_mat+1;
        
        m_z1(i1,i2,:)=sum(reshape(I_1block,[Bsq,Nbu]),1);
        m_z2(i1,i2,:)=sum(reshape(I_2block,[Bsq,Nbu]),1);
        corr_z1z2(i1,i2,:)=sum(reshape(hyp2n3D_scalar(I_1block,I_2block),[Bsq,Nbu]),1);
%         corr_z1z2(i1,i2,:)=sum(reshape(onion_mult2D(I_1block,I_2block),[Bsq,Nbu]),1);
        corr_z1(i1,i2)=sum(reshape(sum(I_1block.^2,3),[Bsq,1]));
        corr_z2(i1,i2)=sum(reshape(sum(I_2block.^2,3),[Bsq,1]));
    end
    fprintf(sprintf('%s%%3d%%%%', repmat('\b', 1, 4)), round(i1/L1f*100));
end
fprintf('. Done!\n')

%Final step to evaluate expected values
m_z1=m_z1./Bsq;
m_z2=m_z2./Bsq;
corr_z1z2=corr_z1z2./Bsq;
corr_z1=corr_z1./Bsq;
corr_z2=corr_z2./Bsq;

%Calculation of Q2^n index
m2_z1=sum(m_z1.^2,3);
m2_z2=sum(m_z2.^2,3);
var_z1=corr_z1-m2_z1;
var_z2=corr_z2-m2_z2;
cov_z1z2=sqrt(sum((corr_z1z2-hyp2n3D_scalar(m_z1,m_z2)).^2,3));
% cov_z1z2=sqrt(sum((corr_z1z2-onion_mult2D(m_z1,m_z2)).^2,3));

mean_bias=2*sqrt(m2_z1.*m2_z2)./(m2_z1+m2_z2);
var_sum=var_z1+var_z2;
index=(var_sum==0);
var_sum(index)=eps;
Q2n_map=2*cov_z1z2./var_sum.*mean_bias;
Q2n_map(index)=mean_bias(index);
Q2n_index=mean2(Q2n_map);

%Map resizing to match original image sizes
%Q2n_map=imresize(Q2n_map,[L1,L2],'nearest');

end