function [ I_out ] = interp_7thorderpp( I_in, b )
%INTERP_SINC_ODD Odd sinc interpolation
%   I_in is the input set of low resolution images
%   b is the upscaling ratio

h=load_pp_filter(b);
Lh=length(h);
E=floor(Lh/2); % edge extension
I_in=edge_extend(I_in,E);
[L1,L2,Nb]=size(I_in);

I_out1=zeros(L1*b+2*floor(Lh/2)-2*floor(b/2),L2,Nb);
rep_h=repmat(reshape(h,[Lh,1,1]),[1,L2,Nb]);
for i1=1:L1
    iA=(i1-1)*b+(1:Lh);
    I_out1(iA,:,:)=I_out1(iA,:,:)+repmat(I_in(i1,:,:),[Lh,1,1]).*rep_h;
end
I_out1=I_out1(E*b+floor(Lh/2)-floor(b/2)+1:end-E*b-floor(Lh/2)+floor(b/2),:,:);
L1b=size(I_out1,1);
I_out=zeros(L1b,L2*b+2*floor(Lh/2)-2*floor(b/2),Nb);
rep_h=repmat(reshape(h,[1,Lh,1]),[L1b,1,Nb]);
for i2=1:L2
    iA=(i2-1)*b+(1:Lh);
    I_out(:,iA,:)=I_out(:,iA,:)+repmat(I_out1(:,i2,:),[1,Lh,1]).*rep_h;
end
I_out=I_out(:,E*b+floor(Lh/2)-floor(b/2)+1:end-E*b-floor(Lh/2)+floor(b/2),:);

