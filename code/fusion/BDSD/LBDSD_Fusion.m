function MS_out = LBDSD_Fusion( Pan, MS_in, ratio, GNyq_PAN, GNyq_MS, B1, B2 )
% LBDSD_FUSION  Local BDSD MMSE Fusion
%   Pan is the panchromatic image
%   MS_in is the multispectral upsampled set of images
%   ratio is the scale ratio between original and upsampled MS image set
%   GNyq_PAN is the gain of the PAN sensor MTF at Nyquist frequency
%   GNyq_MS is the gain of the MS sensor MTF at Nyquist frequency
%   B1 is the height of the block used to segmentate the image
%   B2 is its length

if nargin<4
    GNyq_PAN=0.15;
end
if nargin<5
    GNyq_MS= 0.29 .* ones(1,size(MS_in,3));
end

[L1_FR,L2_FR,Nb]=size(MS_in);
rem_L1=rem(L1_FR,ratio);
rem_L2=rem(L2_FR,ratio);
if rem_L1~=0
    L1_FR=L1_FR-rem(L1_FR,ratio);
    MS_in=MS_in(1:L1_FR,:,:);
    Pan=Pan(1:L1_FR,:);
end
if rem_L2~=0
    L2_FR=L2_FR-rem(L2_FR,ratio);
    MS_in=MS_in(:,1:L2_FR,:);
    Pan=Pan(:,1:L2_FR);
end

L1=L1_FR/ratio;
L2=L2_FR/ratio;
MS_in=double(MS_in);
Pan=double(Pan);

P_L=MTF_filter(Pan,GNyq_PAN,ratio);
P_L=P_L(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end);
MS_in_LR=imresize(MS_in,1/ratio);
MS_in_L=MTF_filter(MS_in_LR,GNyq_MS,ratio);

L1f=ceil(L1/B1); edge_ud=(L1f*B1-L1)/2; edge_ud_FR=edge_ud*ratio;
L2f=ceil(L2/B2); edge_lr=(L2f*B2-L2)/2; edge_lr_FR=edge_lr*ratio;
edge_extension=[floor(edge_ud),floor(edge_lr);ceil(edge_ud),ceil(edge_lr)];
edge_extension_FR=[floor(edge_ud_FR),floor(edge_lr_FR);ceil(edge_ud_FR),ceil(edge_lr_FR)];
P_L=edge_extend(P_L,edge_extension);
MS_in_LR=edge_extend(MS_in_LR,edge_extension);
MS_in_L=edge_extend(MS_in_L,edge_extension);
Pan=edge_extend(Pan,edge_extension_FR);
MS_in=edge_extend(MS_in,edge_extension_FR);

Bsq=B1*B2;
B1_FR=B1*ratio;
B2_FR=B2*ratio;
Bsq_FR=B1_FR*B2_FR;
MS_out=zeros(L1f*B1*ratio,L2f*B2*ratio,Nb);
for i1=1:L1f
    iA=(i1-1)*B1+(1:B1);
    iA_FR=(i1-1)*B1_FR+(1:B1_FR);
    for i2=1:L2f
        iB=(i2-1)*B2+(1:B2);
        iB_FR=(i2-1)*B2_FR+(1:B2_FR);
        
        MS_in_Lblock=reshape(MS_in_L(iA,iB,:),[Bsq,Nb]);
        Hd=cat(2,MS_in_Lblock,reshape(P_L(iA,iB),[Bsq,1]));
        gamma=pinv(Hd)*(reshape(MS_in_LR(iA,iB,:),[Bsq,Nb])-MS_in_Lblock);
        %gamma=(Hd'*Hd)\Hd'*(reshape(MS_in_LR(iA,iB,:),[Bsq,Nb])-MS_in_Lblock);
        
        MS_inblock=reshape(MS_in(iA_FR,iB_FR,:),[Bsq_FR,Nb]);
        H=cat(2,MS_inblock,reshape(Pan(iA_FR,iB_FR),[Bsq_FR,1]));
        MS_out_temp=MS_inblock+H*gamma;
        MS_out(iA_FR,iB_FR,:)=reshape(MS_out_temp,[B1_FR,B2_FR,Nb]);
    end
end

MS_out=MS_out(floor(edge_ud_FR)+(1:L1*ratio),floor(edge_lr_FR)+(1:L2*ratio),:);

end
