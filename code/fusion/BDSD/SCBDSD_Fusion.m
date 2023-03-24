function [MS_out, index_map, K] = SCBDSD_Fusion( Pan, MS_in, scale, GNyq_PAN, GNyq_MS, index_map, K, B )
%BDSD_FUSION  BDSD MMSE Fusion
%   Pan is the panchromatic image
%   MS_in is the multispectral upsampled set of images
%   MS_in_LR is its original resolution set of multispectral images
%   scale is the scale between original and upsampled MS image set

if nargin<4
    GNyq_PAN=0.15;
end
if nargin<5
    GNyq_MS= 0.29 .* ones(1,size(MS_in,3));
end
if nargin==6
    K=max(index_map(:));
elseif nargin<=5
    K=30;
end
if nargin<=7
    B=31;
end

[L1,L2,Nb]=size(MS_in);

MS_in=double(MS_in);
Pan=double(Pan);

P_L=MTF_filter(Pan,GNyq_PAN,scale);
MS_in_L=MTF_filter(MS_in,GNyq_MS,scale);

if nargin<=5
    Sw=std_block(P_L,B*scale);
    Sw=Sw(:);
    P_Lv=P_L(:);
    index_map=kmeans(cat(2,P_Lv/max(P_Lv),Sw/max(Sw)),K);
end

dim=L1*L2;

MS_in_Lblock=reshape(MS_in,[dim,Nb]);
Hd=cat(2,MS_in_Lblock,ones(dim,1),reshape(P_L,[dim,1]));
gamma_global=pinv(Hd)*(reshape(MS_in,[dim,Nb])-MS_in_Lblock);

MS_out=zeros(L1,L2,Nb);
count=0;
for ii=1:K
    index=(index_map==ii);
    index_MS=repmat(index,[1,1,Nb]);
    Nc=nnz(index);
    MS_in_Lblock=reshape(MS_in_L(index_MS),[Nc,Nb]);
    Hd=cat(2,MS_in_Lblock,ones(Nc,1),reshape(P_L(index),[Nc,1]));
    gamma=pinv(Hd)*(reshape(MS_in(index_MS),[Nc,Nb])-MS_in_Lblock);
    if any(gamma(Nb+2,:)<0)
        gamma=gamma_global;
        count=count+1;
    end
    
    MS_inblock=reshape(MS_in(index_MS),[Nc,Nb]);
    H=cat(2,MS_inblock,ones(Nc,1),reshape(Pan(index),[Nc,1]));
    MS_out(index_MS)=MS_inblock+H*gamma;
end
if count~=0
    disp(['Suppressed clusters: ',num2str(count)]);
end

index_map=reshape(index_map,[L1,L2]);
end

