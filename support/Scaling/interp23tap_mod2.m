function I_Interpolated = interp23tap_mod2(I_Interpolated,ratio,order,flag_align,flag_edge)
% keyboard
if nargin<=4
    flag_edge='circular';  % Also available 'symmetric','replicate'
end
if nargin<=3
    flag_align=0;
end
if nargin<=2 || isempty(order)
    order = 11; % Lagrange polynomial interpolation order
end

[r,c,b] = size(I_Interpolated);

BaseCoeff = intfilt(ratio,order,'Lagrange');

I1LRU = zeros(ratio.*r, ratio.*c, b); 
if flag_align==0 || rem(ratio,2)~=0
    I1LRU(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:) = I_Interpolated;
else
    I1LRU(floor(ratio/2):ratio:end,floor(ratio/2):ratio:end,:) = I_Interpolated;
end

for ii = 1 : b
    t = I1LRU(:,:,ii);
    t = imfilter(t',BaseCoeff,flag_edge); 
    I1LRU(:,:,ii) = imfilter(t',BaseCoeff,flag_edge); 
end

I_Interpolated = I1LRU;

end