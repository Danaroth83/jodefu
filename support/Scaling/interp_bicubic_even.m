function [ I_out ] = interp_bicubic_even( I_in, r )
%INTERP_BICUBIC_EVEN Even bicubic interpolation
%   I_in is the input low resolution image
%   r is the upscaling ratio (if odd, there is 1/2 pixel top left disalignment)

h=load_cubic_filter(r,'e');
Lh=length(h);
hsq=h'*h;
E=2; %Put E=0 if no edge treatment is desired
I_in=edge_extend(I_in,E);
[L1,L2,Nb]=size(I_in);
I_out=zeros((L1-1)*r+Lh,(L2-1)*r+Lh,Nb);
for k=1:Nb
    for i1=1:L1
        for i2=1:L2
            I_out((i1-1)*r+1:(i1-1)*r+Lh,(i2-1)*r+1:(i2-1)*...
                r+Lh,k)=I_out((i1-1)*r+1:(i1-1)*r+Lh,(i2-1)*...
                r+1:(i2-1)*r+Lh,k)+I_in(i1,i2,k)*hsq;
        end
    end
end
I_out=I_out(floor((Lh-r)/2+E*r+1):floor(end-(Lh-r)/2-E*r),...
    floor((Lh-r)/2+E*r+1):floor(end-(Lh-r)/2-E*r),:);