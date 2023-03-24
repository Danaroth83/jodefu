function [ h ] = load_pp_filter( ratio )
%LOAD_PP_FILTER Loads piecewise polynomial interpolation filter
% Reference: Image Reconstruction by Convolution with Symmetrical Piecewise
% nth-Order Polynomial Kernels [1999: Meijering, Zuiderveld, Viergever]

len=8;

hlen=len/2;
alpha=-71/83232;
a=zeros(hlen,len);
a(1,:)=[254*alpha+821/1734,-621*alpha-1148/867,0,760*alpha+1960/867,0,-384*alpha-1393/578,0,1];
a(2,:)=[301*alpha+1687/6936,-3309*alpha-2492/867,14952*alpha+32683/2312,-35640*alpha-128695/3468,...
    47880*alpha+127575/2312,-36000*alpha-13006/289,14168*alpha+120407/6936,-2352*alpha-2233/1156];
a(3,:)=[57*alpha+35/6936,-1083*alpha-175/1734,8736*alpha+1995/2312,-38720*alpha-4725/1156,...
    101640*alpha+1575/136,-157632*alpha-5670/289,133336*alpha+42525/2312,-47280*alpha-8505/1156];
a(4,:)=[alpha,-27*alpha,312*alpha,-2000*alpha,7680*alpha,-17664*alpha,22528*alpha,-12288*alpha];

if mod(ratio,2)==0
    x=(0:hlen-1)'*ones(1,ratio)+ones(hlen,1)*(1/(2*ratio):1/ratio:1-1/(2*ratio));
else
    x=(0:hlen-1)'*ones(1,ratio)+ones(hlen,1)*(0:1/ratio:1-1/ratio);
end
h=zeros(1,ratio*hlen);
for ii=1:hlen
    t=repmat(x(ii,:)',[1,len]).^(repmat(len-1:-1:0,[ratio,1]));
    h((ii-1)*ratio+(1:ratio))=sum(repmat(a(ii,:),[ratio,1]).*t,2);
end
if mod(ratio,2)==0
    h=[fliplr(h),h];
else
    h=[fliplr(h(2:end)),h];
end
