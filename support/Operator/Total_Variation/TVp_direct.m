%% New Total Variation Operator
%
% Description:
% This function applies the digital total variation operator described in
% [1], and concatenates the 6 obtained gradiants over the 4th dimension.
% It assumes that the input was obtained from TV_direct.m
% Specifically, given a gradient x=[x1,x2], where x1 is the gradient over
% the columns and x2 over the row, the six concatenated output gradients 
% are defined as [1]:
%     - Column gradient 1: y(:,:,1)=  x1
%     - Column gradient 2: y(m,n,2)= (x1(m,n)+x1(m-1,n)+x1(m,n+1)+x1(m-1,n+1))/4
%     - Column gradient 3: y(m,:,3)= (x1(m,:)+x1(m-1,:))/2
%     - Row gradient 1:    y(:,:,4)=  x2
%     - Row gradient 2:    y(m,n,5)= (x2(m,n)+x2(m+1,n)+x2(m,n-1)+x2(m+1,n-1))/4
%     - Row gradient 3:    y(:,n,6)= (x2(:,n)-x2(:,n-1))/2
%
%
% Usage:
% y=TVd_direct(x);
%
% Input:
% x: An Nd matrix who has the column and row gradients over the 4th dimension
%
% Output:
% y: The calculated 6 gradients concatenated over the 4th dimension
%
% References:
% [1] Condat L.,  "Discrete total variation: New definition and
% minimization." SIAM Journal on Imaging Sciences 10.3 (2017): 1258-1290.

function y=TVp_direct(x)

    % x=TV_direct(x);
    
    sizes=size(x);
    x=reshape(x,sizes(1),sizes(2),sizes(3),sizes(4),[]);

    y=zeros(size(x,1),size(x,2),size(x,3),6,size(x,5));    
    y(:,:,:,1,:)=x(:,:,:,1,:);
    y(2:end,1:end-1,:,2,:)=(x(2:end,1:end-1,:,1,:)+x(1:end-1,1:end-1,:,1,:)+x(2:end,2:end,:,1,:)+x(1:end-1,2:end,:,1,:))/4; 
    y(1,1:end-1,:,2,:)=(x(1,1:end-1,:,1,:)+x(1,2:end,:,1,:))/4; 
    y(2:end,:,:,3,:) = (x(2:end,:,:,1,:)+x(1:end-1,:,:,1,:))/2;
    y(1,:,:,3,:) = x(1,:,:,1,:)/2;
    
    y(:,:,:,4,:)=x(:,:,:,2,:);
    y(1:end-1,2:end,:,5,:)=(x(2:end,1:end-1,:,2,:)+x(1:end-1,1:end-1,:,2,:)+x(2:end,2:end,:,2,:)+x(1:end-1,2:end,:,2,:))/4; 			
    y(1:end-1,1,:,5,:)=(x(1:end-1,1,:,2,:)+x(2:end,1,:,2,:))/4;
    y(:,2:end,:,6,:) = (x(:,2:end,:,2,:)+x(:,1:end-1,:,2,:))/2; 
    y(:,1,:,6,:) = x(:,1,:,2,:)/2;

    y=reshape(y,[sizes(1),sizes(2),sizes(3),6,sizes(5:end)]);
    
end

%% Alternative implementation (slower, also returns some different result in non relevant positions)

% function y = TVp_direct2(x)
%     
%     % x=TV_direct(x);
%     
%     sizes=size(x);
%     x=reshape(x,sizes(1),sizes(2),sizes(3),sizes(4),[]);
%    
%     x=padarray(x,[1,1],0,'both');
%     
%     y=zeros([size(x,1),size(x,2),size(x,3),6,size(x,5)]);
%     y(:,:,:,[1,4],:)=x;
%     y(:,:,:,2,:)=(x(:,:,:,1,:)+circshift(x(:,:,:,1,:),1,1)+circshift(x(:,:,:,1,:),-1,2)+circshift(circshift(x(:,:,:,1,:),-1,2),1,1))/4; 
%     y(:,:,:,5,:)=(x(:,:,:,2,:)+circshift(x(:,:,:,2,:),1,2)+circshift(x(:,:,:,2,:),-1,1)+circshift(circshift(x(:,:,:,2,:),-1,1),1,2))/4; 
%     y(:,:,:,3,:)=(x(:,:,:,1,:)+circshift(x(:,:,:,1,:),1,1))/2;
%     y(:,:,:,6,:)=(x(:,:,:,2,:)+circshift(x(:,:,:,2,:),1,2))/2;
%     
%     y([1,end],:,:,:,:)=[]; y(:,[1,end],:,:,:)=[];
%     
%     y=reshape(y,[sizes(1),sizes(2),sizes(3),6,sizes(5:end)]);
%     
% end