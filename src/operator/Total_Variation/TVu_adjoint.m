%% Upwind Total Variation Adjoint Operator
%
% Description:
% This function applies the adjoint over the upwind total variation [1]
% operation over a matrix in the format of the output of TVu_direct.m
% The four concatenated gradients of the upwind Total Variation are defined
% as:
%     - Column gradient 1: y(n,:,1)=max(x(n,:)-x(n-1,:),0)
%     - Row gradient 1:    y(:,n,2)=max(x(:,n)-x(:,n-1),0)
%     - Column gradient 2: y(n,:,3)=max(x(n+1,:)-x(n,:),0)
%     - Row gradient 2:    y(:,n,4)=max(x(:,n+1)-x(:,n),0)
% Notice the 0 lower thresholding is not implemented in this version of the
% script, to preserve the linearity of the operator
%
% Usage:
% y=TVu_adjoint(x);
%
% Input:
% x: An Nd matrix whose upwind TV gradients are on the 4th dimension
%
% Output:
% y: The result of the adjoint operation
%
% References:
% [1] Chambolle A., Levine S. E. and Lucier B. J., "An upwind 
%     finite-difference method for total variation–based image smoothing.", 
%     SIAM Journal on Imaging Sciences 4.1 (2011): 277-299.

function u=TVu_adjoint(u)

    u=ipermute(u,[2:4,1,5:ndims(u)]);
    sizes=size(u);
    u=reshape(u,sizes(1),sizes(2),sizes(3),[]);
    
    u =   cat(2, u(1,1,:,:),            diff(u(1,:,:,:),1,2))...
         +cat(3, u(2,:,1,:),            diff(u(2,:,:,:),1,3))...
         +cat(2, -diff(u(3,:,:,:),1,2), u(3,end,:,:))...
         +cat(3, -diff(u(4,:,:,:),1,3), u(4,:,end,:));
    
    u=reshape(u,sizes(2:end));
    
end

% %% Original version (Just works on 3D Matrices)
%
% function u=TVu_adjoint(u)
%     u = [u(1,:,:,1);diff(u(:,:,:,1),1,1)]+[u(:,1,:,2) diff(u(:,:,:,2),1,2)]+...
% 		  [-diff(u(:,:,:,3),1,1);u(end,:,:,3)]+[-diff(u(:,:,:,4),1,2) u(:,end,:,4)];	
% end