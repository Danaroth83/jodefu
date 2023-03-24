function output=R2(I_PAN,I_MS_LR,ratio)

if nargin<=2, ratio=size(I_PAN,1)/size(I_MS_LR,1); end

I_MS_LR = double(I_MS_LR);
I_PAN = double(I_PAN);

[L1,L2,Nb]=size(I_MS_LR);
I_PAN_LR=imresize(I_PAN,1/ratio);

IHc = reshape(I_PAN_LR,[L1*L2 1]);
ILRc = reshape(I_MS_LR,[L1*L2 Nb]);
w = ILRc\IHc;

I=sum(I_MS_LR.*repmat(reshape(w,[1,1,Nb]),[L1,L2,1]),3);

output=1-var(I_PAN_LR(:)-I(:))/var(I_PAN_LR(:));

end


