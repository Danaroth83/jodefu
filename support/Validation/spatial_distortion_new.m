function [Ds Ds_map] = spatial_distortion_new(fusa,pan,ref_ridotta,pan_ridotta,m2,ratio,Dim)
% [Ds Ds_map] = spatial_distortion_new(fusa,pan,ref_ridotta,pan_ridotta,m2,ratio,Dim);
%
% keyboard
fusa = double(fusa);
pan = double(pan);
ref_ridotta = double(ref_ridotta);
pan_ridotta = double(pan_ridotta);

N3 = size(fusa,3);

aux1=vett_spat_mappa(fusa,pan,Dim);
aux2=vett_spat_mappa(ref_ridotta,pan_ridotta,Dim/ratio);

d=aux1-aux2;

Ds = ((1/N3)*sum(abs(d).^m2,3)).^(1/m2);
Ds_map = imresize(Ds, Dim, 'nearest');
% figure('Name','Spatial Distortion','NumberTitle','off'), imshow(Ds_map,[]),title('Spatial Distortion')
Ds = mean2(Ds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aux=vett_spat_mappa(dato,pan,Dim);

%%%CONTROLLO NEL CASO DI IMMAGINI QUADRATE E NON MULTIPLE DI DIM%%%%
[N1,N2,N3]=size(dato);
stepx=ceil(N1/Dim);
stepy=ceil(N2/Dim);

if stepy<=0
    stepy=1;
    stepx=1;
end;

est1=stepx*Dim-N1;
est2=stepy*Dim-N2;

if sum([(est1~=0),(est2~=0)])>0
  datoo=[];
  for k=1:N3
  appo=zeros(N1+est1,N2+est2);
  appo(1:N1,1:N2)=squeeze(dato(:,:,k));
  appo(:,N2+1:N2+est2)=appo(:,N2:-1:N2-est2+1);
  appo(N1+1:N1+est1,:)=appo(N1:-1:N1-est1+1,:);
  datoo=cat(3,datoo,appo);
  end
  dato=datoo;
  clear datoo
  %[N1,N2,N3]=size(dato)
  
  appo=zeros(N1+est1,N2+est2);
  appo(1:N1,1:N2)=pan;
  appo(:,N2+1:N2+est2)=appo(:,N2:-1:N2-est2+1);
  appo(N1+1:N1+est1,:)=appo(N1:-1:N1-est1+1,:);
  pan=appo;
  clear appo
  
end

N3=size(dato,3);
aux = zeros(stepx,stepy,N3);
for i=1:N3
    aux(:,:,i) = img_qinew_mappa(squeeze(dato(:,:,i)),pan,Dim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mappa=img_qinew_mappa(a,b,sizer)

shift=sizer;
[N1,N2]=size(a);
stepx=ceil(N1/shift);
stepy=ceil(N2/shift);

if stepy<=0
    stepy=1;
    stepx=1;
end;
valori=zeros(stepx,stepy);

for j=1:stepx
    for i=1:stepy
        rit_ima=a(((j-1)*shift)+1:((j-1)*shift)+sizer,((i-1)*shift)+1:((i-1)*shift)+sizer);
        rit_imb=b(((j-1)*shift)+1:((j-1)*shift)+sizer,((i-1)*shift)+1:((i-1)*shift)+sizer);
   
        valori(j,i)=img_qi(rit_ima,rit_imb,sizer);
    end;
end; 
mappa=valori;

