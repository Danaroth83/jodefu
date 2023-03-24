function I1LR = interp23tap(I1LR,Resize_fact)

if (2^round(log2(Resize_fact)) ~= Resize_fact)
    disp('Error: Only resize factors power of 2');
    return;
end 

[r,c,b] = size(I1LR);

CDF23 = 2.*[0.5 0.305334091185 0 -0.072698593239 0 0.021809577942 0 -0.005192756653 0 0.000807762146 0 -0.000060081482];
CDF23 = [fliplr(CDF23(2:end)) CDF23];
BaseCoeff = CDF23;
first = 1;

for z = 1 : Resize_fact/2

    I1LRU = zeros((2^z) * r, (2^z) * c, b);
    
    if first
        I1LRU(2:2:end,2:2:end,:) = I1LR;
        first = 0;
    else
        I1LRU(1:2:end,1:2:end,:) = I1LR;
    end

    for ii = 1 : b
        t = I1LRU(:,:,ii); 
        t = imfilter(t',BaseCoeff,'circular'); 
        I1LRU(:,:,ii) = imfilter(t',BaseCoeff,'circular'); 
    end
    
    I1LR = I1LRU;
    
end

end