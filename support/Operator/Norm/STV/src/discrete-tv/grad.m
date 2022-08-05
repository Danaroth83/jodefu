function [gx,gy] = grad(u)
%% compute (gx,gy) = discrete gradient of u
    gx = u(:,[2:end,end])-u; 
    gy = u([2:end,end],:)-u; 
end

