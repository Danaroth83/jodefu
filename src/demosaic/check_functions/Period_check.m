function [mask,I_in,period]=Period_check(mask,I_in)

    %% Test for Periodic Mask
    tol=10E-9*max(mask(:));
    mask_idxnonzero=mask>tol;
    sum_mask_idxnonzero=sum(mask_idxnonzero,3);

    if any(sum_mask_idxnonzero(:)>=2) || any(sum_mask_idxnonzero(:)<=0)
        error('Only non overlapping and fully covering masks are allowed');
    end
    
    I_in=I_in./sum(mask,3);
    mask=double(mask_idxnonzero);
    period=period_custom(mask);
    
    if isempty(period)
        error('The mask is not periodic');
    end
    
%     mask_new=zeros(size(mask,1),size(mask,2),prod(period));
%     perm=zeros(1,prod(period));
%     for ii=1:period(1)
%         for jj=1:period(2)
%             idx=(ii-1)*period(2)+jj;
%             mask_new(ii:period(1):end,jj:period(2):end,idx)=1;
%             perm(idx)=find(mask(ii,jj,:)==1);
%         end
%     end
    
%    options.perm=perm;
%    options.flip=[0,0];
%    options.edge=[0,0];

end