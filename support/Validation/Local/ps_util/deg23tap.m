function If = deg23tap(If,Resize_fact)

if (2^round(log2(Resize_fact)) ~= Resize_fact)
    disp('Error: Only resize factors power of 2');
    return;
end 

[~,~,b] = size(If);

CDF23 = [0.5 0.305334091185 0 -0.072698593239 0 0.021809577942 0 -0.005192756653 0 0.000807762146 0 -0.000060081482];
CDF23 = [fliplr(CDF23(2:end)) CDF23];
BaseCoeff = CDF23;
first = 1;

for z = 1 : Resize_fact/2

    for ii = 1 : b
        t = If(:,:,ii); 
        t = imfilter(t',BaseCoeff,'replicate'); 
        h(:,:,ii) = imfilter(t',BaseCoeff,'replicate'); 
    end
    
    if first
        h = h(2:2:end,2:2:end,:);
        first = 0;
    else
        h = h(1:2:end,1:2:end,:);
    end

    If = h;
end

end