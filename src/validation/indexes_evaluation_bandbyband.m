function [SCC_index,Q_index,ERGAS_index]=indexes_evaluation_bandbyband(I_1,I_2,ratio)
Nb=size(I_1,3);
SCC_index=zeros(1,Nb);
Q_index=zeros(1,Nb);
ERGAS_index=zeros(1,Nb);

for ii=1:Nb
    SCC_index(ii)=SCC(I_1(:,:,ii),I_2(:,:,ii),'laplacian','global1');
    Q_index(ii)=Q(I_1(:,:,ii),I_2(:,:,ii));
    ERGAS_index(ii)=ERGAS(I_1(:,:,ii),I_2(:,:,ii),ratio);
end

end