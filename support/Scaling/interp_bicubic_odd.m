function [ I_out ] = interp_bicubic_odd( I_in, b )
%INTERP_BICUBIC_ODD Odd bicubic interpolation
%   I_in is the input set of low resolution images
%   b is the upscaling ratio (if even, there is 1/2 pixel top left disalignment)

h=load_cubic_filter(b,'o');
Lh=length(h);
hsq=h'*h;
E=2; %Put E=0 if no edge treatment is desired
I_in=edge_extend(I_in,E);
[L1,L2,Nb]=size(I_in);
I_out=zeros((L1-1)*b+Lh,(L2-1)*b+Lh,Nb);
for k=1:Nb
    for i1=1:L1
        for i2=1:L2
            I_out((i1-1)*b+1:(i1-1)*b+Lh,(i2-1)*b+1:(i2-1)*...
                b+Lh,k)=I_out((i1-1)*b+1:(i1-1)*b+Lh,(i2-1)*...
                b+1:(i2-1)*b+Lh,k)+I_in(i1,i2,k)*hsq;
        end
    end
end

I_out=I_out(floor((Lh-b)/2+E*b+1):floor(end-(Lh-b)/2-E*b),...
    floor((Lh-b)/2+E*b+1):floor(end-(Lh-b)/2-E*b),:);