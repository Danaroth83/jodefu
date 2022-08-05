function [ sCC,sCC_block ] = SCC( I_1, I_2 , preproc, type, Q_blocks, K )
%SCC Calculates the spatial cross correlation between I_1 and I_2
% I_1 and I_2 must be same-sized image sets
% preproc: edge enhancement procedure ('sobel2','sobel','laplacian','prewitt2','prewitt','none') 
% type identifies if the SCC is global or with moving blocks
%   ('global1' removes the outmost pixel after edge enhanchement)
% K are the overlapping pixels of moving blocks

if nargin<=2 || isempty(preproc), preproc='laplacian'; end
if nargin<=3 || isempty(type), type='global1'; end
if nargin<=4 || isempty(Q_blocks), Q_blocks=32; end
if nargin<=5 || isempty(K), K=1; end

if strcmp(preproc,'none')
    I_1p=I_1;
    I_2p=I_2;
    if strcmpi(type,'global1'), type='global'; end
elseif strcmp(preproc,'sobel2')|| strcmp(preproc,'prewitt2')
    preproc2=preproc(1:end-1);
    fpreproc2=fspecial(preproc2);
    I_1p = zeros(size(I_1));
    for idim=1:size(I_1,3)
        % I_1p(:,:,idim)= imfilter(I_1(:,:,idim),fpreproc2);
        % I_1p(:,:,idim)= imfilter(I_1p(:,:,idim),fpreproc2');
        I_1p(:,:,idim)= sqrt(imfilter(I_1(:,:,idim),fpreproc2).^2+imfilter(I_1(:,:,idim),fpreproc2').^2);
    end
    I_2p = zeros(size(I_2));
    for idim=1:size(I_2,3)
        % I_2p(:,:,idim)= imfilter(I_2(:,:,idim),fpreproc2);
        % I_2p(:,:,idim)= imfilter(I_2p(:,:,idim),fpreproc2');
        I_2p(:,:,idim)=sqrt(imfilter(I_2(:,:,idim),fpreproc2).^2+imfilter(I_2(:,:,idim),fpreproc2').^2);
    end
else
    fpreproc=fspecial(preproc);
    I_1p = zeros(size(I_1));
    for idim=1:size(I_1,3)
        I_1p(:,:,idim)= imfilter(I_1(:,:,idim),fpreproc);
    end
    I_2p = zeros(size(I_2));
    for idim=1:size(I_2,3)
        I_2p(:,:,idim)= imfilter(I_2(:,:,idim),fpreproc);
    end
end

if strcmpi(type,'global1') && min(size(I_1,1),size(I_1,2))>=3
    I_1p=I_1p(2:end-1,2:end-1,:);
    I_2p=I_2p(2:end-1,2:end-1,:);
end


if strncmpi(type,'global',6)
    I_1_v=I_1p(:);
    I_2_v=I_2p(:);
    sCC=sum(I_1_v.*I_2_v);
    sCC=sCC/sqrt(sum(I_1_v.^2));
    sCC=sCC/sqrt(sum(I_2_v.^2));
else
    edge_extension=[floor((Q_blocks-K)/2);ceil((Q_blocks-K)/2)];
    I_1p=edge_extend(I_1p,edge_extension);
    I_2p=edge_extend(I_2p,edge_extension);
    [L1e,L2e,Nb]=size(I_1p);
    L1f=floor((L1e-Q_blocks)/K+1);
    L2f=floor((L2e-Q_blocks)/K+1);
    Bsq=Q_blocks^2;
    I1sq=I_1p.^2;
    I2sq=I_2p.^2;
    I1I2=I_1p.*I_2p;
    m1_block=zeros(L1f,L2f);
    m2_block=zeros(L1f,L2f);
    sp_block=zeros(L1f,L2f);
    for i1=1:L1f
        iA=(i1-1)*K+(1:Q_blocks);
        for i2=1:L2f
            iB=(i2-1)*K+(1:Q_blocks);
            m1_block(i1,i2)=sum(reshape(I1sq(iA,iB,:),[Bsq*Nb,1]));
            m2_block(i1,i2)=sum(reshape(I2sq(iA,iB,:),[Bsq*Nb,1]));
            sp_block(i1,i2)=sum(reshape(I1I2(iA,iB,:),[Bsq*Nb,1]));
        end
    end
    sCC_block=sp_block./sqrt(m1_block)./sqrt(m2_block);
    sCC=sum(sCC_block(:))/numel(sCC_block);
end

end

