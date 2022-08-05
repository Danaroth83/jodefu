function [ sCC,sCC_map ] = SCC_mask( I1, I2 ,preproc, type, mask )
%SCC Calculates the spatial cross correlation between I_1 and I_2
% I_1 and I_2 must be same-sized image sets
% preproc: edge enhancement procedure ('sobel2','sobel','laplacian','prewitt2','prewitt','none') 
% type identifies if the SCC is global or with moving blocks
%   ('global1' removes the outmost pixel after edge enhanchement)
% K are the overlapping pixels of moving blocks

if nargin<=2 || isempty(preproc), preproc='laplacian'; end
if nargin<=3 || isempty(type), type='global'; end
if nargin<=4 || isempty(mask), mask=ones(size(I1,1),size(I1,2)); end

if strcmp(preproc,'none')
    I1_p=I1;
    I2_p=I2;
elseif strcmp(preproc,'sobel2')|| strcmp(preproc,'prewitt2')
    preproc2=preproc(1:end-1);
    I1_p = zeros(size(I1));
    for idim=1:size(I1,3),
        I1_p(:,:,idim)= imfilter(I1(:,:,idim),fspecial(preproc2));
        I1_p(:,:,idim)= imfilter(I1_p(:,:,idim),(fspecial(preproc2))');
    end
    I2_p = zeros(size(I2));
    for idim=1:size(I2,3),
        I2_p(:,:,idim)= imfilter(I2(:,:,idim),fspecial(preproc2));
        I2_p(:,:,idim)= imfilter(I2_p(:,:,idim),(fspecial(preproc2))');
    end
else
    I1_p = zeros(size(I1));
    for idim=1:size(I1,3),
        I1_p(:,:,idim)= imfilter(I1(:,:,idim),fspecial(preproc));
    end
    I2_p = zeros(size(I2));
    for idim=1:size(I2,3),
        I2_p(:,:,idim)= imfilter(I2(:,:,idim),fspecial(preproc));
    end
end

if strcmpi(type,'global1')
    I1_p=I1_p(2:end-1,2:end-1,:);
    I2_p=I2_p(2:end-1,2:end-1,:);
    mask=mask(2:end-1,2:end-1);
end

Nc=max(mask(:));
sCC=zeros(1,Nc);
sCC_map=zeros(size(I1_p,1),size(I1_p,2));

for nc=1:Nc
    mask_temp=(mask==nc);
    Nm=nnz(mask_temp);

    I1_v=zeros(Nm,size(I1_p,3));
    I2_v=zeros(Nm,size(I2_p,3));
    for ii=1:size(I1_p,3)
        I1_temp=I1_p(:,:,ii);
        I2_temp=I2_p(:,:,ii);
        I1_v(:,ii)=I1_temp(mask_temp);
        I2_v(:,ii)=I2_temp(mask_temp);
    end

    sCC_map(mask_temp)=sum(I1_v.*I2_v,2)*Nm/sqrt(sum(I1_v(:).^2))/sqrt(sum(I2_v(:).^2));
    sCC(nc)=sum(sCC_map(mask_temp))/Nm;
end

sCC_map(mask<=0)=NaN;

end

