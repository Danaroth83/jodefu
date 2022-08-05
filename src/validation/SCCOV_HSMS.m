function sCC = SCCOV_HSMS( I_1, I_2MS, flag_histmatch )
%SCC Calculates the spatial cross correlation between I_1 and I_2
% I_1 is HS
% I_2 is MS
if nargin<=2
    flag_histmatch=0;
end

for ii=1:size(I_1,3)
    I_1(:,:,ii)=I_1(:,:,ii)-mean2(I_1(:,:,ii));
end
for ii=1:size(I_2MS,3)
    I_2MS(:,:,ii)=I_2MS(:,:,ii)-mean2(I_2MS(:,:,ii));
end

[L1,L2,Nb]=size(I_1);
Nb_MS=size(I_2MS,3);
I_1_rep=repmat(reshape(I_1,[L1*L2,Nb,1]),[1,1,Nb_MS]);
if flag_histmatch==0
    I_2_rep=repmat(reshape(I_2MS,[L1*L2,1,Nb_MS]),[1,Nb,1]);
else
    I_2_rep=reshape(I_2MS,[L1*L2,Nb,Nb_MS]);
end
sCC=squeeze(sum(I_1_rep.*I_2_rep,1)./sqrt(sum(I_1_rep.^2,1).*sum(I_2_rep.^2,1)));

% if strcmp(type,'global')
%     I_1_v=I_1(:);
%     I_2_v=I_2(:);
% else
%     edge_extension=[floor((Q_blocks-K)/2);ceil((Q_blocks-K)/2)];
%     I_1=edge_extend(I_1,edge_extension);
%     I_2=edge_extend(I_2,edge_extension);
%     [L1e,L2e,Nb]=size(I_1);
%     L1f=floor((L1e-Q_blocks)/K+1);
%     L2f=floor((L2e-Q_blocks)/K+1);
%     Bsq=Q_blocks^2;
%     I1sq=I_1.^2;
%     I2sq=I_2.^2;
%     I1I2=I_1.*I_2;
%     m1_block=zeros(L1f,L2f);
%     m2_block=zeros(L1f,L2f);
%     sp_block=zeros(L1f,L2f);
%     for i1=1:L1f
%         iA=(i1-1)*K+(1:Q_blocks);
%         for i2=1:L2f
%             iB=(i2-1)*K+(1:Q_blocks);
%             m1_block(i1,i2)=sum(reshape(I1sq(iA,iB,:),[Bsq*Nb,1]));
%             m2_block(i1,i2)=sum(reshape(I2sq(iA,iB,:),[Bsq*Nb,1]));
%             sp_block(i1,i2)=sum(reshape(I1I2(iA,iB,:),[Bsq*Nb,1]));
%         end
%     end
%     sCC_block=sp_block./sqrt(m1_block)./sqrt(m2_block);
%     sCC=sum(sCC_block(:))/numel(sCC_block);
% end

end

