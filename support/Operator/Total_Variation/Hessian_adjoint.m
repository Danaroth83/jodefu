%% Hessian Adjoint Operator
%
% Description:
% This function calculates the adjoint of Hessian operator given the 
% 4 gradients calculated by Hessian_direct.m, and returns the column and 
% row gradient over the 4th dimension
% The Hessian operator is describes as follows: 
% given x=[x1,x2], where x1 is the gradient over the columns and x2 over 
% the row, the four concatenated output gradients are defined as [1]:
%     - Column Hessian:     y(n,:,1)=x1(n,:)-x1(n-1,:)
%     - Column/Row Hessian: y(:,n,2)=x1(:,n)-x1(:,n-1)
%     - Row Hessian:        y(:,n,3)=x2(:,n)-x2(:,n-1)
%     - Row/Column Hessian: y(n,:,4)=x2(n,:)-x2(n-1,:)
%
% Usage:
% y=Hessian_adjoint(x);
%
% Input:
% x: An Nd matrix whose 4 Hessians are concateneted over the
%    4th dimension
%
% Output:
% y: The result of the adjoint operator
%
% References:
% [1] Bredies K., Kunisch K., and Pock T., "Total generalized variation." 
% SIAM Journal on Imaging Sciences 3.3 (2010): 492-526.

function u=Hessian_adjoint(u)
    
    u=ipermute(u,[2:4,1,5:ndims(u)]);
    sizes=size(u);
    u=reshape(u,sizes(1),sizes(2),sizes(3),[]);
    
    u= cat(1, cat(2,-diff(u(1,:,:,:),1,2),u(1,end,:,:)) - cat(3,u(2,:,1,:),diff(u(2,:,:,:),1,3)),...
		      cat(3,-diff(u(3,:,:,:),1,3),u(3,:,end,:)) - cat(2,u(4,1,:,:),diff(u(4,:,:,:),1,2)));
          
    u=reshape(u,[2,sizes(2:end)]);
    u=permute(u,[2:4,1,5:ndims(u)]);
    
    % u=TV_adjoint(u);

end

% %% Original implementation (works with just 3D matrices)
%
% function u=Hessian_adjoint(u)
%     u= cat(4,...
% 		 [-diff(u(:,:,:,1),1,1);u(end,:,:,1)]-[u(:,1,:,2) diff(u(:,:,:,2),1,2)],...
% 		 [-diff(u(:,:,:,3),1,2) u(:,end,:,3)]-[u(1,:,:,4);diff(u(:,:,:,4),1,1)]);
% end