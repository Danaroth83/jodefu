%% BINARY MASK CREATION
%  
% Author:
% Daniele Picone
%
% Description:
% This script creates a mask according to the binary mask creation
% suggested in [1] (no algorithmic procedure is described in the 
% paper, but this script implements the rationale described in it), 
% with minimum periodicity.
%
% Usage:
% mask=binarytree_mask(Nb,mask_type)
%
% Input:
% Nb: Amount of bands for the mask
% mask_type: Method for the generation of the mask, among the following
%       'u': Assigns samples with the best available uniformity
%       'd': Assigns a half pixels to a dominant band, and uniformly the rest
%       's': Assigns pixels with hierarchic dominance (First band has half
%            pixels, second band has 1/4 and so on)
%
% Output: 
% mask: A 2D matrix filled with the indices associated to the
%         assigned band
%
% Reference:
% [1] Monno Y., Kikuchi S., Tanaka M. and Okutomi M., "A Practical One-Shot
%     Multispectral Imaging System Using a Single Image Sensor", IEEE TIP
%     2015 Oct; 24(10):3048-59.

function mask_MS=binarytree_mask(Nb,masktype)
    if Nb==1, mask_MS=1; return; end
    mask=1;
    if strncmpi(masktype,'u',1)     % Uniform binary mask
        cycles=nextpow2(Nb);
        for ii=1:cycles
            mask=binarytree_recursion(mask,Nb);
        end
    elseif strncmpi(masktype,'d',1) % Dominant band binary mask
        cycles=nextpow2(Nb-1);
        mask=1;
        mask=binarytree_recursion(mask,2);
        mask1=mask(:,:,1);
        for ii=1:cycles
            mask1=binarytree_recursion(mask1,Nb-1);
        end
        mask=cat(3,mask1,repmat(mask(:,:,2),[size(mask1,1)/size(mask,1),size(mask1,2)/size(mask,2),1]));
    elseif strncmpi(masktype,'s',1) % Sequence of dominant bands
        mask=1;
        for ii=1:Nb-1
            mask_old=mask(:,:,2:end);
            mask=binarytree_recursion(mask(:,:,1),2);
            mask=cat(3,mask,repmat(mask_old,[size(mask,1)/size(mask_old,1),size(mask,2)/size(mask_old,2),1]));
        end
    end
    mask_MS=zeros(size(mask,1),size(mask,2));
    for kk=1:size(mask,3), mask_MS(mask(:,:,kk)==1)=size(mask,3)-kk+1; end
end



function mask_out=binarytree_recursion(mask,N)
    mask_out=[];
    for kk=1:min(N-size(mask,3),size(mask,3))
        mask_current=mask(:,:,kk);
        if sum(mask_current(:))==1
            mask1=[zeros(size(mask_current)),mask_current;mask_current,zeros(size(mask_current))];
            mask2=[mask_current,zeros(size(mask_current));zeros(size(mask_current)),mask_current];
            mask_out=cat(3,mask_out,mask1,mask2);
        else
            idx_ones=find(mask_current(:)==1);
            mask1=zeros(numel(mask_current),1);
            mask1(idx_ones(2))=1;
            mask2=zeros(numel(mask_current),1);
            mask2(idx_ones(1))=1;
            mask_out=cat(3,mask_out,reshape(mask1,size(mask_current)),reshape(mask2,size(mask_current)));
        end
    end
    mask=repmat(mask,[size(mask_out,1)/size(mask,1),size(mask_out,2)/size(mask,2),1]);
    mask_out=cat(3,mask_out,mask(:,:,min(N-size(mask,3))+1:end));
end