%% Proximal operator of the l_2,infinity,1 norm
% 
% Description:
% This function calculates the proximal operator of the l_{2,inf,1} norm, 
% defined as:
%   l2inf1norm = @(x) gamma * sum(max(sqrt(sum(x.^2,1)),[],2),3:ndims(x));
% where gamma>=0
% This is equivalent to calculating y=x-projball_l21inf(x,gamma), where:
%   projball_l21inf(x,gamma) = argmin_{x' : l21infnorm(x')<=gamma} ||x_x'||_2
% is the projection of of x over the l_2,1,inf ball of radius gamma.
% Also notice that the proximal operator for the l_{inf,2,1} norm can be
% obtained as:
% prox_linf21(x,gamma)=ipermute(prox_l2inf1(permute(x,[2,1,3:ndims(x)]),gamma),[2,1,3:ndims(x)]);
%
% Usage:
% y = prox_2inf1(x, gamma);
%
% Input:
% x: Input Matrix, whose norm l_2 is applied on the first dimension,
%    l_inf on the the second,l_1 on the rest
% gamma: scalar multiplier for the norm
% 
% Output:
% y: Evaluation of the proximal operator

function y = prox_l2inf1(x, gamma)
    tmp1 = sqrt(sum(x.^2,1));
    tmp2 = min(max(max((cumsum(sort(tmp1,2,'descend'),2)-gamma)./(1:size(x,2)),[],2),0)./tmp1,1);
    y = x .* tmp2;
end

% %% Alternative implementation
% function y = prox_l2inf1 (x, gamma)
% 	
%     sizex=size(x);
%     x=reshape(x,size(x,1),size(x,2),[]);
%     
%     y=zeros(size(x));
%     tmp = sqrt(sum(x.^2,1));
%     mask = sum(tmp,2)>gamma;
%     tau = max((cumsum(sort(tmp,2,'descend'))-gamma)./(1:size(x,2)),[],2);
%     tmp=tmp(:,:,mask);
% 	  y(:,:,mask) = x(:,:,mask) .* (tau./max(tmp, tau));
%     
%     y=reshape(y,sizex);
%     
% end


% %% Alternative version (only works for gamma>0)
% function y = prox_l2inf1(x, gamma)
%     tmp1=sqrt(sum(x.^2,1));
%     v=projball_l1inf(shiftdim(tmp1,1),gamma);
%     y=x.*max(tmp1-shiftdim(v,-1),0)./tmp1;
%     y(isnan(y))=0;
% end