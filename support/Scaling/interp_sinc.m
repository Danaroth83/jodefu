function [ I_out ] = interp_sinc( I_in, b )
%INTERP_SINC_ODD Odd sinc interpolation
%   I_in is the input set of low resolution images
%   b is the upscaling ratio (if even, there is 1/2 pixel top left disalignment)

h=load_sinc_filter(b);
% h=load_cubic_filter(b);
Lh=length(h);
E=floor(Lh/2); % edge extension
I_in=edge_extend(I_in,E);
[L1,L2,Nb]=size(I_in);

%% method 1
% hsq=h'*h;
% I_out=zeros(L1*b+2*floor(Lh/2)-2*floor(b/2),L2*b+2*floor(Lh/2)-2*floor(b/2),Nb);
% rep_hsq=repmat(hsq,[1,1,Nb]);
% for i1=1:L1
%     iA=(i1-1)*b+(1:Lh);
%     for i2=1:L2
%         iB=(i2-1)*b+(1:Lh);
%         I_out(iA,iB,:)=I_out(iA,iB,:)+repmat(I_in(i1,i2,:),[Lh,Lh,1]).*rep_hsq;
%     end
% end
% I_out=I_out(E*b+floor(Lh/2)-floor(b/2)+1:end-E*b-floor(Lh/2)+floor(b/2),E*b+floor(Lh/2)-floor(b/2)+1:end-E*b-floor(Lh/2)+floor(b/2),:);

%% method 2
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