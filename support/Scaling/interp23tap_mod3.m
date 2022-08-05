function I_Interpolated = interp23tap_mod3(I_Interpolated,ratio,order,flag_align,flag_edge)
% keyboard
if numel(ratio)==1, ratio=[ratio,ratio]; end
if nargin<=4, flag_edge='circular'; end % Also available 'symmetric','replicate'
if nargin<=3, flag_align=0; end
if nargin<=2 || isempty(order)
    order = 11; % Lagrange polynomial interpolation order
end

[r,c,b] = size(I_Interpolated);

BaseCoeff_x = intfilt(ratio(1),order,'Lagrange');
BaseCoeff_y = intfilt(ratio(2),order,'Lagrange');

I1LRU = zeros(ratio(1).*r, ratio(2).*c, b);
if flag_align==0 || rem(ratio(1),2)~=0,
    idx_x=floor(ratio(1)/2)+1:ratio(1):size(I1LRU,1);
else
    idx_x=floor(ratio(1)/2):ratio(1):size(I1LRU,1);
end
if flag_align==0 || rem(ratio(2),2)~=0,
    idx_y=floor(ratio(2)/2)+1:ratio(2):size(I1LRU,2);
else
    idx_y=floor(ratio(2)/2):ratio(2):size(I1LRU,2);
end
I1LRU(idx_x,idx_y,:) = I_Interpolated;

for ii = 1 : b
    t = I1LRU(:,:,ii);
    t = imfilter(t',BaseCoeff_x,flag_edge); 
    I1LRU(:,:,ii) = imfilter(t',BaseCoeff_y,flag_edge); 
end

I_Interpolated = I1LRU;

end