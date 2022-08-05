function [MS_out, index_map_FR, index_map, K] = CBDSD_Fusion( Pan, MS_in, ratio, GNyq_PAN, GNyq_MS, index_map_FR, index_map, K, B )
%BDSD_FUSION  BDSD MMSE Fusion
%   index_map_FR is the cluster map at high resolution
%   Pan is the panchromatic image
%   MS_in is the multispectral upsampled set of images
%   ratio is the scale ratio between original and upsampled MS image set
%   GNyq_PAN is the gain of the PAN sensor MTF at Nyquist frequency
%   GNyq_MS is the gain of the MS sensor MTF at Nyquist frequency
%   K is the number of clusters
%   B is the block size for calculating local standard deviation statistics

if nargin<=3
    GNyq_PAN=0.15;
end
if nargin<=4
    GNyq_MS= 0.29 .* ones(1,size(MS_in,3));
end
if nargin==6
    index_map=imresize(index_map_FR,1/ratio,'nearest');
end
if nargin==6 || nargin==7
    K=max(index_map_FR(:));
elseif nargin<=5
    K=30;
end
if nargin<=8
    B=31;
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

P_L_FR=MTF_filter(Pan,GNyq_PAN,ratio);
P_L=P_L_FR(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end);
MS_in_LR=imresize(MS_in,1/ratio);
MS_in_L=MTF_filter(MS_in_LR,GNyq_MS,ratio);

dim=L1*L2;

if nargin<=5
    Sw=std_block(P_L,B);
    Sw=Sw(:);
    P_Lv=P_L(:);
    
    % index_map=reshape(kmeans(cat(2,P_Lv/max(P_Lv),Sw/max(Sw)),K),[L1,L2]);
    % index_map_FR=imresize(index_map,scale,'nearest');
    [index_map,C]=kmeans(cat(2,P_Lv/max(P_Lv),Sw/max(Sw)),K);
    %index_map_FR=imresize(index_map,scale,'nearest');

    index_map_FR=zeros(L1_FR,L2_FR);
    scale_sq=ratio^2;
    for ii=1:scale_sq
        i1=rem(ii-1,ratio)+1;
        i2=(ii-i1)/ratio+1;
        % [i1,i2]=ind2sub(ii,[scale,scale]);
        iA=i1:ratio:L1_FR;
        iB=i2:ratio:L2_FR;
        P_Lcut=P_L_FR(iA,iB);
        Sw=std_block(P_Lcut,B);
        P_Lcutv=P_Lcut(:);
        Sw=Sw(:);
        v_in=cat(2,P_Lcutv/max(P_Lcutv),Sw/max(Sw));
        % index_temp=kmeans(v_in,K,'start',C,'MaxIter',0);
        [~,index_temp]=min(squeeze(sum((repmat(reshape(v_in',[1,2,dim]),[K,1,1])-repmat(C,[1,1,dim])).^2,2)));
        % for i3=1:dim
        %     [~,index_temp(ii)]=min(sum((repmat(v_in(i3,:),[K,1])-C).^2,2));
        % end
        index_map_FR(iA,iB)=reshape(index_temp,[L1,L2]);
    end
end

MS_in_Lblock=reshape(MS_in_L,[dim,Nb]);
Hd=cat(2,MS_in_Lblock,ones(dim,1),reshape(P_L,[dim,1]));
gamma_global=pinv(Hd)*(reshape(MS_in_LR,[dim,Nb])-MS_in_Lblock);

MS_out=zeros(L1*ratio,L2*ratio,Nb);
count=0;
for ii=1:K
    index=(index_map==ii);
    index_MS=repmat(index,[1,1,Nb]);
    Nc=nnz(index);
    index_FR=(index_map_FR==ii);
    index_MS_FR=repmat(index_FR,[1,1,Nb]);
    Nc_FR=nnz(index_FR);
    
    MS_in_Lblock=reshape(MS_in_L(index_MS),[Nc,Nb]);
    Hd=cat(2,MS_in_Lblock,ones(Nc,1),reshape(P_L(index),[Nc,1]));
    gamma=pinv(Hd)*(reshape(MS_in_LR(index_MS),[Nc,Nb])-MS_in_Lblock);
    if any(gamma(Nb+2,:)<0)
        gamma=gamma_global;
        count=count+1;
    end
    
    MS_inblock=reshape(MS_in(index_MS_FR),[Nc_FR,Nb]);
    H=cat(2,MS_inblock,ones(Nc_FR,1),reshape(Pan(index_FR),[Nc_FR,1]));
    MS_out(index_MS_FR)=MS_inblock+H*gamma;
end
if count~=0
    disp(['Suppressed clusters: ',num2str(count)]);
end
end

