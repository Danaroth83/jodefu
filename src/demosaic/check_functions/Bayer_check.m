function [mask,I_in,options]=Bayer_check(mask,I_in)

    flip=[0,0];
    edge=[0,0];
    %% Test for Bayer Mask
    tol=10E-9*max(mask(:));
    mask_idxnonzero=mask>tol;
    sum_mask_idxnonzero=sum(mask_idxnonzero,3);

    if any(sum_mask_idxnonzero(:)>=2) || any(sum_mask_idxnonzero(:)<=0)
        error('Only non overlapping and fully covering masks are allowed');
    end
    
    I_in=I_in./sum(mask,3);
    mask=double(mask_idxnonzero);
    if ~isequal(mask(1:end-2,1:end-2,:),mask(3:end,3:end,:))
        error('Mask hasn''t the proper periodicity');
    end
    if size(mask,3)>=4 || size(mask,3)<=1
        error('Mask isn''t in Bayer form');
    end
    
    mask_period=mask(1:2,1:2,:);
    gband=[];
    for kk=1:size(mask,3)
        if mask_period(1,1,kk)==1 && mask_period(2,2,kk)==1
            gband=kk;
            break;
        elseif mask_period(1,2,kk)==1 && mask_period(2,1,kk)==1
            gband=kk;
            if mod(size(mask,1),2)==0
                flip=[1,0];
            elseif mod(size(mask,2),2)==0
                flip=[0,1];
            else
                flip=[1,0]; edge=[1,0];
            end
            break;
        end       
    end
    if isempty(gband), error('The Mask is not Bayer applicable'); end
    
    if edge(1)==1
        mask=mask([1:size(mask,1),1],:,:);
        I_in=I_in([1:size(mask,1),1],:,:); 
    end
    if flip(1)==1
        mask=mask(end:-1:1,:,:);
        I_in=I_in(end:-1:1,:,:); 
    end
    if flip(2)==1
        mask=mask(:,end:-1:1,:);
        I_in=I_in(:,end:-1:1,:); 
    end
    
    mask_period=mask(1:2,1:2,:);
    for kk=1:size(mask,3)
        if mask_period(1,2,kk)==1, rband=kk; end
        if mask_period(2,1,kk)==1, bband=kk; end
    end
    
    perm=[rband,gband,bband];
    % I_in=I_in(:,:,perm);
    mask=mask(:,:,perm);
    
    options.perm=perm;
    options.flip=flip;
    options.edge=edge;

end