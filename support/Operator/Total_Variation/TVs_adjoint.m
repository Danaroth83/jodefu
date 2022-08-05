% ADJOINT FOR THE SHANNON TOTAL VARIATION
%
% Wrapper of sdiv
function y=TVs_adjoint(x,n)
%     sizes=size(x);
%     x=reshape(x,size(x,1),size(x,2),size(x,3),size(x,4),[]);
%     y=zeros(size(x,1)/n,size(x,2)/n,size(x,3),size(x,5));
%     for ii=1:size(y,4)
%         y(:,:,:,ii)=sdiv(x(:,:,:,:,ii),n);
%     end
%     y=reshape(y,[size(y,1),size(y,2),size(y,3),sizes(5:end)]);
      
    x=ipermute(x,[1:3,ndims(x),4:ndims(x)-1]); 
    sizes=size(x);
    x=reshape(x,sizes(1),sizes(2),[],2);
    y=zeros(size(x,1)/n,size(x,2)/n,size(x,3));
    x=x/n^2; %% is this correct?
    for ii=1:size(x,3)
        y(:,:,ii)=-sdiv(x(:,:,ii,1),x(:,:,ii,2),n); %% Compared to Moisan code, minus is necessary to define a proper adjoint operator
    end
    y=reshape(y,[size(y,1),size(y,2),sizes(3:end-1)]);

end