function [map,gradient] = grad(im,method)

if isequal(method,'supremum')==1
    [gradient,map] = supgrad(im);
elseif isequal(method,'metric')==1
    [gradient,map] = mbgrad(im);
elseif isequal(method,'robust')==1
    [gradient,map] = rcmgrad(im);
else
    error('incorrect input method, choose between supremum, metric or robust')
end