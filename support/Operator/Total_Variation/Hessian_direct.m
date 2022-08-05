%% Hessian Operator
%
% Description:
% This function calculates the Hessian operator given the gradients of the 
% first two dimension of a Nd matrix (it supposes to have the result of 
% TV_direct.m as input). 
% Specifically, given a gradient x=[x1,x2], where x1 is the gradient over
% the columns and x2 over the row, the four concatenated output gradients 
% are defined as [1]:
%     - Column Hessian:     y(n,:,1)=x1(n,:)-x1(n-1,:)
%     - Column/Row Hessian: y(:,n,2)=x1(:,n)-x1(:,n-1)
%     - Row Hessian:        y(:,n,3)=x2(:,n)-x2(:,n-1)
%     - Row/Column Hessian: y(n,:,4)=x2(n,:)-x2(n-1,:)
%
% Usage:
% y=TV_direct(x);
%
% Input:
% x: An Nd matrix whose row and column gradients are concateneted over the
%    4th dimension
%
% Output:
% y: The four calculated gradients concatenated over the 4th dimension
%
% References:
% [1] Bredies K., Kunisch K., and Pock T., "Total generalized variation." 
% SIAM Journal on Imaging Sciences 3.3 (2010): 492-526.

function x=Hessian_direct(x)
    
    % x=TV_direct(x);
    
    x=ipermute(x,[2:4,1,5:ndims(x)]);
    sizes=size(x); 
    x=reshape(x,sizes(1),sizes(2),sizes(3),[]);
    
    x = cat(1, cat(2, x(1,1,:,:)          , diff(x(1,:,:,:),1,2)          ),...
               cat(3, diff(x(1,:,:,:),1,3), zeros(1,sizes(2),1,size(x,4)) ),...
               cat(3, x(2,:,1,:)          , diff(x(2,:,:,:),1,3)          ),...
               cat(2, diff(x(2,:,:,:),1,2), zeros(1,1,sizes(3),size(x,4)) )); 

    x=reshape(x,[4,sizes(2:end)]);
    x=permute(x,[2:4,1,5:ndims(x)]);
    
end

%% Original implementation (just works with 3d matrices)
%
% function x=Hessian_direct(x)
%     x= cat(3,...
% 		    [x(1,:,:,1);diff(x(:,:,:,1),1,1)],[diff(x(:,:,:,1),1,2) zeros(size(x,1),1,size(x,3))],...
% 		    [x(:,1,:,2) diff(x(:,:,:,2),1,2)],[diff(x(:,:,:,2),1,1);zeros(1,size(x,2),size(x,3))]);
% end