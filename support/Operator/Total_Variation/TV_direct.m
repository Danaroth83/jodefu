%% Total Variation Operator
%
% Description:
% This function applies the total variation operator over a color image.
% Specifically it calculates the gradient over the first two dimension and
% concatenates the result over along the 4th dimension of a given Nd matrix
% The two concatenated gradients are defined as:
%     - Column gradient: y(n,:)=x(n,:)-x(n-1,:)
%     - Row gradient:    y(:,n)=x(:,n)-x(:,n-1)
%
% Usage:
% y=TV_direct(x);
%
% Input:
% x: An Nd matrix with at least two dimensions
%
% Output:
% y: The given gradients concatenated over the 4th dimension

function x=TV_direct(x)
    
    sizes=size(x);
    x=shiftdim(x,-1);

    x=cat(1,cat(2, diff(x,1,2), zeros([1,1,sizes(2:end)])         ),...
            cat(3, diff(x,1,3), zeros([1,sizes(1),1,sizes(3:end)])));

    x=permute(x,[2:4,1,5:ndims(x)]);

end

% %% Original version (Just works on 2D or 3D Matrices)
%
% function x=TV_direct(x)
%     x=cat(4,[diff(x,1,1);zeros(1,size(x,2),size(x,3))],[diff(x,1,2) zeros(size(x,1),1,size(x,3))]);
% end