%% CREATION OF COMPOSED HR AND LR IMAGE MASK
%
% Description:
% The function creates a mask composed mask of a High Resolution (HR) and
% Low Resolution (LR) mask, by inserting the LR mask values in the zero
% valued samples of the HR mask
%
% Usage:
% mask=create_combo_mask(HR_mask,LR_mask,N1,N2);
%
% Example:
% HR_mask=[1,1;1,0]; LR_mask=[1,2;3,4];
% mask=create_combo_mask(HR_mask,LR_mask);
%
% Input:
% HR_mask: A binary mask, the zero values have to be filled with LR mask
%          values
% LR_mask: A mask whose value indicates which LR band the final mask will
%          be assigned to (eg: 3 indicates it is assigned to the third band)
%          (Note: values marked with zero will be assigned randomly,
%                 values marked with -1 will be have a random combination of input with sum 1)
% N1: Bands of the HR image (default: max(HR_mask(:)))
% N2: Bands of the LR image (default: max(LR_mask(:)))
%
% Output:
% mask: A binary mask, whose bands are N1+N2, the first N1 bands are
%       assigned to the HR image, and the last N2 to the LR one; the
%       spatial sides are automatically chosen to take into account the
%       periodicity of the combined contributions

function mask=create_combo_mask(pan_mask,ms_mask,sizes_final_mask,Nbp)

%% Parsing inputs

if nargin<=0, pan_mask=0; end
if nargin<=1, ms_mask=[1,2;3,4]; end
if nargin<=2, sizes_final_mask=[]; end
if nargin<=3, Nbp=max(pan_mask(:)); end

if isempty(sizes_final_mask), Nbm=max(ms_mask(:)); else, Nbm=sizes_final_mask(3); end


%% Initialization of the masks

pan_mask_idx=pan_mask;
pan_mask_idx(pan_mask~=0)=1; % Find all non-zero positions in the HR mask

Nbp2=ceil(max(pan_mask(:)));
ms_mask(ms_mask>0)=Nbp2+ms_mask(ms_mask>0);

%% Combination of the mask

if numel(ms_mask)==0        % Case of empty mask
    mask=[]; return;
elseif numel(ms_mask)==1    % Case of LR mask with a single pixel
    pan_mask(pan_mask==0)=ms_mask;
else 
    % Selection of all rows and columns fully covered by the HR mask
    sum_col=size(pan_mask,1)-sum(pan_mask_idx,1);
    sum_row=size(pan_mask,2)-sum(pan_mask_idx,2);
    col=(sum_col==0);
    row=(sum_row==0);
    
    % Finding the dimension of the mask without the pixels assigned to the HR mask
    col_L1=max(sum_col(~col));
    if col_L1==0, warning('No LR signal is included in the mask'); end
    row_L2=max(sum_row(~row));
    col_L2=sum(~col);
    row_L1=sum(~row); 
    
    cond_col=all(sum_col(:,~col)==col_L1);
    cond_row=all(sum_row(~row)==row_L2);
    
    % Case 1: number of pixels to assign over the column is constant
    if (cond_col && cond_row && row_L1>col_L1) || (cond_col && ~cond_row) 
        mul_L1=lcm(col_L1,size(ms_mask,1));
        mul_L2=lcm(col_L2,size(ms_mask,2));
        pan_mask=repmat(pan_mask,[mul_L1/col_L1,mul_L2/col_L2]);
        ms_mask=repmat(ms_mask,[mul_L1/size(ms_mask,1),mul_L2/size(ms_mask,2)]);
        
        pan_mask(pan_mask==0)=ms_mask(:);
     
    % Case 2: number of pixels to assign over the row is constant
    elseif (cond_col && cond_row && row_L1<=col_L1) || (~cond_col && cond_row)
        mul_L1=lcm(row_L1,size(ms_mask,1));
        mul_L2=lcm(row_L2,size(ms_mask,2));
        pan_mask=repmat(pan_mask,[mul_L1/row_L1,mul_L2/row_L2]);
        ms_mask=repmat(ms_mask,[mul_L1/size(ms_mask,1),mul_L2/size(ms_mask,2)]);
        
        pan_mask=pan_mask.';
        ms_mask=ms_mask.';
        pan_mask(pan_mask==0)=ms_mask(:);
        pan_mask=pan_mask.';
        
    % Case 3: Assigning the empty pixels without keeping the structure
    else
        inputval=sum(pan_mask(:)==0);
        ms_mask=repmat(ms_mask(:),ceil(inputval/numel(ms_mask)));
        pan_mask(pan_mask==0)=ms_mask(1:inputval);
    end
end

%% Extending the mask to the desired sizes
if ~isempty(sizes_final_mask)
    %pan_mask=wextend(2,'ppd',pan_mask,[sizes_final_mask(1)-size(pan_mask,1),sizes_final_mask(2)-size(pan_mask,2)],'rr');
    pan_mask=repmat(pan_mask,[ceil(sizes_final_mask(1)/size(pan_mask,1)),ceil(sizes_final_mask(2)/size(pan_mask,2))]);
    pan_mask=pan_mask(1:sizes_final_mask(1),1:sizes_final_mask(2),:);
end

mask=zeros([size(pan_mask,1),size(pan_mask,2),Nbp+Nbm]);

%% Assigning the pixel randomly when ms_mask=0;
num_zeros=sum(pan_mask(:)==0);
idx_zeros=find(pan_mask(:)==0);
sample_length=[round(num_zeros/Nbm)*ones(1,Nbm-1),num_zeros-round(num_zeros/Nbm)*(Nbm-1)];
for jj=1:Nbm
    idx_ms=randsample(length(idx_zeros),sample_length(jj));
    pan_mask(idx_zeros(idx_ms))=jj+Nbp2;
    idx_zeros(idx_ms)=[];
end

mask=reshape(mask,[],Nbp+Nbm);

%% Assigning the binary valued pixels when ms_mask=-1;
num_minusone=sum(pan_mask(:)==-1);
sample_length=round(num_minusone/2);
for jj=Nbp+(1:Nbm)
    mask(randsample(num_minusone,sample_length),jj)=1/(Nbp+Nbm);
end

%% Assigning the randomly weighted pixels when ms_mask=-2;
mask=reshape(mask,[],Nbp+Nbm);
ind_minustwo=find(pan_mask(:)==-2);
mask(ind_minustwo,Nbp+(1:Nbm))=simplex_uniform(Nbm,length(ind_minustwo));

mask=reshape(mask,size(pan_mask,1),size(pan_mask,2),[]);

%% Generating the mask from the positions assigned

for jj=1:min(Nbp2,Nbp)
    mask(:,:,jj)=mask(:,:,jj)+(pan_mask==jj);
end

for jj=Nbp2+(1:Nbm)
    mask(:,:,jj+Nbp-Nbp2)=mask(:,:,jj+Nbp-Nbp2)+(pan_mask==jj);
end


