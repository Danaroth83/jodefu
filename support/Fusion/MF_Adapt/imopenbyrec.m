function Y = imopenbyrec(I, S)
% Perform opening by reconstruction
%  I: input image
%  S: structuring element (instance of class strel)
%  Y: output image
% TODO: add connectivity as input and do some checks

% S(ceil(size(S,1)/2),ceil(size(S,2)/2)) = 1;

Y = imreconstruct(imerode(I,S),I);