function MS_out = GBDSD_Fusion( Pan, MS_in, ratio, GNyq_PAN, GNyq_MS )
%GBDSD_FUSION  Global BDSD MMSE Fusion
%   Pan is the panchromatic image
%   MS_in is the multispectral upsampled set of images
%   ratio is the scale ratio between original and upsampled MS image set
%   GNyq_PAN is the gain of the PAN sensor MTF at Nyquist frequency
%   GNyq_MS is the gain of the MS sensor MTF at Nyquist frequency

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

% Method 1: taking the original as Low Resolution option (add it as input)
% MS_in_LR=double(MS_in_LR);

% Method 2: Resizing with a close-to-ideal filter for LR option
MS_in_LR=imresize(MS_in,1/ratio);
MS_in_L=MTF_filter(MS_in_LR,GNyq_MS,ratio);

dim=L1*L2;
MS_in_Lblock=reshape(MS_in_L,[dim,Nb]);
Hd=cat(2,MS_in_Lblock,reshape(P_L,[dim,1]));
gamma=pinv(Hd)*(reshape(MS_in_LR,[dim,Nb])-MS_in_Lblock);

dim_FR=dim*ratio^2;
MS_inblock=reshape(MS_in,[dim_FR,Nb]);
H=cat(2,MS_inblock,Pan(:));
MS_out=reshape(MS_inblock+H*gamma,[L1*ratio,L2*ratio,Nb]);

end