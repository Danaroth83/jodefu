%% Total Variation Adjoint Operator
%
% Description:
% This function implements the adjoint operator for the gradient over the
% first two dimension
% It assumes as input a ND matrix, whose gradients are over the 4th
% dimension (The input image has the format of the output of TV_direct)
% The two gradient are defined as:
%     - Column gradient:   y(n,:,1)=x(n,:)-x(n-1,:)
%     - Row gradient:      y(:,n,2)=x(:,n)-x(:,n-1)
%
% Usage:
% y=TV_direct(x);
%
% Input:
% x: An Nd matrix with at least 4 dimension; the 4th dimension (required 
%    size: 2) must contain the vertical and horizontal gradients
%
% Output:
% y: The result of the adjoint operation


function u=TV_adjoint(u)
    
    u=ipermute(u,[2:4,1,5:ndims(u)]);
    sizes=size(u);
    u=reshape(u,sizes(1),sizes(2),sizes(3),[]);
    
    u= -cat(2, u(1,1,:,:), diff(u(1,:,:,:),1,2))...
       -cat(3, u(2,:,1,:), diff(u(2,:,:,:),1,3));
    
    u=reshape(u,sizes(2:end));
    
end

% %% Original implementation (Only valid for 3D Matrices)
% 
% function u=TV_adjoint(u)
%     u=-[u(1,:,:,1);diff(u(:,:,:,1),1,1)]-[u(:,1,:,2) diff(u(:,:,:,2),1,2)];
% end