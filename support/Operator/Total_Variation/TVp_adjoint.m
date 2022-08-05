%% New Total Variation Adjoint Operator
%
% Description:
% This function applies the adjoint of the total variation operator 
% described in [1], and takes the 6 obtained gradients concatenated over
% the 4th dimension to return 2.
% Specifically, given a gradient x=[x1,x2], where x1 is the gradient over
% the columns and x2 over the row, the six concatenated output gradients of
% this total variation are defined as [1]:
%     - Column gradient 1: y(:,:,1)=  x1
%     - Column gradient 2: y(m,n,2)= (x1(m,n)+x1(m-1,n)+x1(m,n+1)+x1(m-1,n+1))/4
%     - Column gradient 3: y(m,:,3)= (x1(m,:)+x1(m-1,:))/2
%     - Row gradient 1:    y(:,:,4)=  x2
%     - Row gradient 2:    y(m,n,5)= (x2(m,n)+x2(m+1,n)+x2(m,n-1)+x2(m+1,n-1))/4
%     - Row gradient 3:    y(:,n,6)= (x2(:,n)-x2(:,n-1))/2
%
%
% Usage:
% y=TVp_adjoint(x);
%
% Input:
% x: An Nd matrix who has the 6 gradients over the 4th dimension
%
% Output:
% y: The result of the adjoint operation
%
% References:
% [1] Condat L.,  "Discrete total variation: New definition and
% minimization." SIAM Journal on Imaging Sciences 10.3 (2017): 1258-1290.

function u=TVp_adjoint(t)
	
    sizes=size(t);
	t=reshape(t,sizes(1),sizes(2),sizes(3),sizes(4),[]);   

    u=zeros([sizes(1:3),2,size(t,5)]);	
    u(1:end-1,2:end,:,1,:)=t(1:end-1,2:end,:,1,:)+(t(1:end-1,2:end,:,2,:)+t(1:end-1,1:end-1,:,2,:)+...
	                         t(2:end,2:end,:,2,:)+t(2:end,1:end-1,:,2,:))/4+(t(1:end-1,2:end,:,3,:)+t(2:end,2:end,:,3,:))/2;
	u(1:end-1,1,:,1,:)=t(1:end-1,1,:,1,:)+(t(1:end-1,1,:,2,:)+t(2:end,1,:,2,:))/4+...
	                   (t(1:end-1,1,:,3,:)+t(2:end,1,:,3,:))/2;
	u(2:end,1:end-1,:,2,:)=t(2:end,1:end-1,:,4,:)+(t(2:end,1:end-1,:,5,:)+t(1:end-1,1:end-1,:,5,:)+...
	                       t(2:end,2:end,:,5,:)+t(1:end-1,2:end,:,5,:))/4+(t(2:end,1:end-1,:,6,:)+t(2:end,2:end,:,6,:))/2;
	u(1,1:end-1,:,2,:)=t(1,1:end-1,:,4,:)+(t(1,1:end-1,:,5,:)+t(1,2:end,:,5,:))/4+...
	                   (t(1,1:end-1,:,6,:)+t(1,2:end,:,6,:))/2;

    u=reshape(u,sizes(1),sizes(2),sizes(3),2,sizes(5:end));
    
    % u=TV_adjoint(u);
end


%% Alternative implementation (slower, also has some nonzero results in non relevant positions)

% function y=TVp_adjoint(x)
% 	
%     sizes=size(x);
% 	  x=reshape(x,sizes(1),sizes(2),sizes(3),sizes(4),[]);
%     
%     x=padarray(x,[1 1],0,'both');
%       
%     y(:,:,:,1,:)=x(:,:,:,1,:)+...
%       (x(:,:,:,2,:)+circshift(x(:,:,:,2,:),1,2)+circshift(circshift(x(:,:,:,2,:),1,2),-1,1)+circshift(x(:,:,:,2,:),-1,1))/4+...
%       (x(:,:,:,3,:)+circshift(x(:,:,:,3,:),-1,1))/2;
%   
%     y(:,:,:,2,:)=x(:,:,:,4,:)+...
%                 (x(:,:,:,5,:)+circshift(x(:,:,:,5,:),1,1)+circshift(circshift(x(:,:,:,5,:),1,1),-1,2)+circshift(x(:,:,:,5,:),-1,2))/4+...
%                 (x(:,:,:,6,:)+circshift(x(:,:,:,6,:),-1,2))/2;
%             
%     y([1,end],:,:,:,:)=[]; y(:,[1,end],:,:,:)=[];
%             
%     y=reshape(y,sizes(1),sizes(2),sizes(3),2,sizes(5:end));
%     
%     % y=TV_adjoint(y);
%                 
% end