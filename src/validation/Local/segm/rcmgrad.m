function [grad,map] = rcmgrad(in)

grad = zeros(size(in(:,:,1),1),size(in(:,:,1),2));
E = padarray(in, [1 1 0], 'symmetric');
[m n] = size(E(:,:,1));

for i=2:m-1
    for j=2:n-1
        N2 = E(i-1:i+1,j-1:j+1,:);
        N2 = reshape(N2,size(N2,3),size(N2,1)^2);
        nn = size(N2,2);
        di = zeros(1,nn*(nn-1)/2);
        cmpt = 1;
        for k=1:nn-1
            x = N2(:,k);
            for l=k+1:nn
                y = N2(:,l);
                di(cmpt) = (sum((x-y).^2)).^0.5;
                cmpt = cmpt+1;
            end
        end
        di = sort(di,'descend');
        grad(i-1,j-1) = di(2);
    end
end

map = double(watershed(grad));