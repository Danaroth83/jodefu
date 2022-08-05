function [varargout] = showNirRG(I, stretch)

if nargin<2
    stretch = [0.01 0.99];
end

if size(I,3) == 8 % worldview 2
    ch = [3,5,8];
else
    ch = [2,3,4];
end

I = viewimage(I(:,:,ch), stretch);

if nargout>0
    varargout{1} = I(:,:,3:-1:1);
end