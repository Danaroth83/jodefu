function [ h ] = load_sinc_filter( b, flag, hlen )
%LOAD_SINC_FILTER Calculates the sinc approximation for FIR interpolation
%filter
%   b is the interpolation scale
%   flag can either be 'e' for even filter or 'o' for odd interpolation
%   hlen is the half length of the filter (full lenght is 2*hlen for even
%       and 2*hlen+1 for odd filter)
if nargin<=1
    if rem(b,2)==0
        flag='e';
    else
        flag='o';
    end
end
if nargin<=2
    hlen=b*11;
end
if flag=='e'
    t=1/(2*b):1/b:1/(2*b)+hlen/b;
    t=[-fliplr(t),t];
    h=sinc(t);
else
   t=1/b:1/b:hlen/b;
   t=[-fliplr(t),0,t];
   h=sinc(t).*kaiser(length(t))';
end

normalize=zeros(1,b);
for ii=1:b
   normalize(ii)=sum(h(ii:b:end));
end
normalize=repmat(normalize,[1,ceil(length(h)/b)]);
h=h./normalize(1:length(h));