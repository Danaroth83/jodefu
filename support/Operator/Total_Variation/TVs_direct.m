%% DIRECT OPERATOR FOR SHANNON TOTAL VARIATION
%
% Wrapper of sgrad
function y=TVs_direct(x,n)
%     sizes=size(x);
%     x=reshape(x,size(x,1),size(x,2),size(x,3),[]);
%     y=zeros([size(x,1)*n,size(x,2)*n,size(x,3),2,size(x,4)]);
%     for ii=1:size(x,4)
%         y(:,:,:,:,ii)=sgrad(x(:,:,:,ii),n);
%     end
%     y=reshape(y,[size(y,1),size(y,2),size(y,3),2,sizes(4:end)]);

      sizes=size(x); if length(sizes)==2, sizes(3)=1; end
      x=reshape(x,size(x,1),size(x,2),[]);
      y=zeros([size(x,1)*n,size(x,2)*n,size(x,3),2]);
      for ii=1:size(x,3)
          [y(:,:,ii,1),y(:,:,ii,2)]=sgrad(x(:,:,ii),n);
      end
      y=y/n^2;
      y=reshape(y,[size(y,1),size(y,2),sizes(3:end),2]);
      y=permute(y,[1:3,ndims(y),4:ndims(y)-1]);
end