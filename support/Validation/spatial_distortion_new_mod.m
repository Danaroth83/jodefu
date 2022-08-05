function [Ds,Q_FR,Q_RR] = spatial_distortion_new_mod(fusa,pan,ref_ridotta,pan_ridotta,m2,Dim,S)
% [Ds Ds_map] = spatial_distortion_new(fusa,pan,ref_ridotta,pan_ridotta,m2,ratio,Dim);
%
% keyboard
if nargin<=6 || isempty(S)
    S=Dim;
end
fusa = double(fusa);
pan = double(pan);
ref_ridotta = double(ref_ridotta);
pan_ridotta = double(pan_ridotta);

N3 = size(fusa,3);

aux1=vett_spat_mappa(fusa,pan,Dim,S);
aux2=vett_spat_mappa(ref_ridotta,pan_ridotta,Dim,S);
[stepx,stepy,~]=size(aux1);
stepsq=stepx*stepy;
Q_FR=sum(reshape(aux1,[stepsq,N3]),1)/stepsq;
Q_RR=sum(reshape(aux2,[stepsq,N3]),1)/stepsq;
Ds = ((1/N3)*sum(abs(Q_FR-Q_RR).^m2)).^(1/m2);
if nargin>=2
    Q_FR=sum(Q_FR)/N3;
    Q_RR=sum(Q_RR)/N3;
end

% d=aux1-aux2;
% Ds = ((1/N3)*sum(abs(d).^m2,3)).^(1/m2);
% Ds_map = imresize(Ds, Dim, 'nearest');
% figure('Name','Spatial Distortion','NumberTitle','off'), imshow(Ds_map,[]),title('Spatial Distortion')
% Ds = mean2(Ds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aux=vett_spat_mappa(I_1,I_2,Dim,S)

edge_extension=[floor((Dim-S)/2);ceil((Dim-S)/2)];
I_1=edge_extend(I_1,edge_extension);
I_2=edge_extend(I_2,edge_extension);
[L1e,L2e,Nb]=size(I_1);
stepx=floor((L1e-Dim)/S+1);
stepy=floor((L2e-Dim)/S+1);

Dimsq=Dim^2;
m1=zeros(stepx,stepy,Nb);
m2=zeros(stepx,stepy);
s1=zeros(stepx,stepy,Nb);
s2=zeros(stepx,stepy);
s12=zeros(stepx,stepy,Nb);
for x=1:stepx
    iA=((x-1)*S)+(1:Dim);
    for y=1:stepy
        iB=((y-1)*S)+(1:Dim);
        I_1block=reshape(I_1(iA,iB,:),[Dimsq,Nb]);
        I_2block=reshape(I_2(iA,iB),[Dimsq,1]);
        m1(x,y,:)=sum(I_1block,1);
        m2(x,y)=sum(I_2block);
        s1(x,y,:)=sum(I_1block.^2,1);
        s2(x,y)=sum(I_2block.^2);
        s12(x,y,:)=sum(I_1block.*repmat(I_2block,[1,Nb]),1);
    end
end
m2r=repmat(m2,[1,1,Nb]);
s2r=repmat(s2,[1,1,Nb]);
mp=m1.*m2r;
msq=m1.^2+m2r.^2;
aux=4*(Dimsq*s12-mp).*mp./(Dimsq*(s1+s2r)-msq)./msq;

