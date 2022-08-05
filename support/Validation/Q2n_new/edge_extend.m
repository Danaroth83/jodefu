function [ I_out ] = edge_extend( I_in, pix_extend, method )
% EDGE_EXTEND Flips the image around the edges
%   I_in is a set of images
%   pix_extend is the numbers of pixel to add; format is [up,left;down,right]
%   As default, if pix_extend is scalar, the extension is considered on every side
%   If it is a 2-variables horizontal array, the extension is considered as [vertical,horizontal]
%   If it is a 2-variables vertical array, the extension is considered as [top&left; bottom&right]
%   method 'r' means the pixels closest to the edge are repeated

if nargin<=2 || method~='r', method='n'; end
if method=='r', K=0; else K=1; end
if numel(pix_extend)==1; pix_extend=[pix_extend,pix_extend; pix_extend,pix_extend]; end
if isequal(size(pix_extend),[1,2]), pix_extend=[pix_extend; pix_extend]; end
if isequal(size(pix_extend),[2,1]), pix_extend=[pix_extend, pix_extend]; end
pixu=pix_extend(1,1); pixl=pix_extend(1,2); pixd=pix_extend(2,1); pixr=pix_extend(2,2);

% Negative extension (it removes extra edges)
if size(I_in,1)<=-(pixu+pixd), I_in=[]; pixu=0; pixd=0; pixl=0; pixr=0; end
if size(I_in,2)<=-(pixl+pixr), I_in=[]; pixu=0; pixd=0; pixl=0; pixr=0; end
if pixl<0, I_in=I_in(:,-pixl+1:end,:); pixl=0; end
if pixr<0, I_in=I_in(:,1:end+pixr,:); pixr=0; end
if pixu<0, I_in=I_in(-pixu+1:end,:,:); pixu=0; end
if pixd<0, I_in=I_in(1:end+pixd,:,:); pixd=0; end

% Recursive procedure if extension is bigger than image sizes
[L1, L2, Nd]=size(I_in);
if L1-K<pixu
    I_in=edge_extend(I_in,[L1-K,0;0,0],method);
    pixu=pixu-L1+K;
    L1=size(I_in,1);
end
if L1-K<pixd
    I_in=edge_extend(I_in,[0,0;L1-K,0],method);
    pixd=pixd-L1+K;
    L1=size(I_in,1);
end
if L2-K<pixl
    I_in=edge_extend(I_in,[0,L2-K;0,0],method);
    pixl=pixl-L2+K;
    L2=size(I_in,2);
end
if L2-K<pixr
    I_in=edge_extend(I_in,[0,0;0,L2-K],method);
    pixr=pixr-L2+K;
    L2=size(I_in,2);
end

I_out=zeros(L1+pixu+pixd,L2+pixl+pixr,Nd);

I_out(pixu+1:end-pixd,pixl+1:end-pixr,:)=I_in;
I_out([1:pixu,end-pixd+1:end],pixl+1:end-pixr,:)=...
    flipud(I_in([end-pixd+1-K:end-K,1+K:pixu+K],:,:));
I_out(pixu+1:end-pixd,[1:pixl, end-pixr+1:end],:)=...
    fliplr(I_in(:,[end-pixr+1-K:end-K,1+K:pixl+K],:));
I_out([1:pixu,end-pixd+1:end],[1:pixl,end-pixr+1:end],:)=rot90(I_in(...
    [end-pixd+1-K:end-K,1+K:pixu+K],[end-pixr+1-K:end-K,1+K:pixl+K],:),2);

end