%% Upwind Total Variation Operator
%
% Description:
% This function applies the upwind total variation [1] operator over a the 
% first two dimension of an Nd matrix.
% Specifically it calculates the gradient over the first two dimension and
% concatenates the result over along the 4th dimension of a given Nd matrix
% The four concatenated gradients are defined as:
%     - Column gradient 1: y(n,:,1)=max(x(n,:)-x(n-1,:),0)
%     - Row gradient 1:    y(:,n,2)=max(x(:,n)-x(:,n-1),0)
%     - Column gradient 2: y(n,:,3)=max(x(n+1,:)-x(n,:),0)
%     - Row gradient 2:    y(:,n,4)=max(x(:,n+1)-x(:,n),0)
% Notice the 0 lower thresholding is not implemented in this version of the
% script, to preserve the linearity of the operator
%
% Usage:
% y=TVu_direct(x);
%
% Input:
% x: An Nd matrix with at least two dimensions
%
% Output:
% y: The given gradients concatenated over the 4th dimension
%
% References:
% [1] Chambolle A., Levine S. E. and Lucier B. J., "An upwind 
%     finite-difference method for total variation–based image smoothing.", 
%     SIAM Journal on Imaging Sciences 4.1 (2011): 277-299.

function x=TVu_direct(x)

    sizes=size(x);
    x=shiftdim(x,-1);
    
    diffx12=diff(x,1,2);
    diffx13=diff(x,1,3);
    zeros3=zeros([1,1,sizes(2:end)]);
    zeros2=zeros([1,sizes(1),1,sizes(3:end)]);
    
    x=cat(1,cat(2, -diffx12, zeros3) ,...
            cat(3, -diffx13, zeros2) ,...
            cat(2, zeros3,   diffx12),...
            cat(3, zeros2,   diffx13));

    x=permute(x,[2:4,1,5:ndims(x)]);
    
end

% %% Original version (Just works on 2D or 3D Matrices)
%
% function x=TVu_direct(x)
%     x=cat(4,[-diff(x,1,1);zeros(1,size(x,2),size(x,3))],[-diff(x,1,2) zeros(size(x,1),1,size(x,3))],...
% 		      [zeros(1,size(x,2),size(x,3));diff(x,1,1)],[zeros(size(x,1),1,size(x,3)) diff(x,1,2)]);
% end