function [ h ] = load_cubic_filter( b, flag )
%LOAD_CUBIC_FILTER Calculates the sinc approximation for FIR interpolation
%filter
%   b is the interpolation scale
%   flag can either be 'e' for even filter or 'o' for odd interpolation
if nargin<=1
    if rem(b,2)==0
        flag='e';
    else
        flag='o';
    end
end
if flag=='e'
    t=1/(2*b):1/b:1-1/(2*b);
    h1=3/2*t.^3-5/2*t.^2+1;
    t=1+1/(2*b):1/b:2-1/(2*b);
    h2=-1/2*t.^3+5/2*t.^2-4*t+2;
    h=[fliplr(h2),fliplr(h1),h1,h2];
else
    t=1/b:1/b:1-1/b;
    h1=3/2*t.^3-5/2*t.^2+1;
    t=1:1/b:2-1/b;
    h2=-1/2*t.^3+5/2*t.^2-4*t+2;
    h=[fliplr(h2),fliplr(h1),1,h1,h2]; 
end
