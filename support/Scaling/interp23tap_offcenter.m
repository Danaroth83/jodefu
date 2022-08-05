function I_in = interp23tap_offcenter(I_in,ratio,order,flag_edge,flag_alignright)
% keyboard
if numel(ratio)==1, ratio=[ratio,ratio]; end
if nargin<=2 || isempty(order), order = 11; end % Lagrange polynomial interpolation order
if nargin<=3, flag_edge='circular'; end         % Also available 'symmetric','replicate'
if nargin<=4, flag_alignright=0; end
if numel(flag_alignright)==1, flag_alignright=[flag_alignright,flag_alignright]; end
if numel(ratio)==1, ratio=[ratio,ratio]; end

BaseCoeff_x=intfilt(ratio(1),order,'Lagrange');
BaseCoeff_y=intfilt(ratio(2),order,'Lagrange');

[r,c,b]=size(I_in);

I1LRU = zeros(ratio(1).*r, ratio(2).*c, b);
if flag_alignright(1)==1 || rem(ratio(1),2)~=0
    idx_x=floor(ratio(1)/2)+1:ratio(1):size(I1LRU,1);
else
    idx_x=floor(ratio(1)/2):ratio(1):size(I1LRU,1);
end
if flag_alignright(2)==1 || rem(ratio(2),2)~=0
    idx_y=floor(ratio(2)/2)+1:ratio(2):size(I1LRU,2);
else
    idx_y=floor(ratio(2)/2):ratio(2):size(I1LRU,2);
end
I1LRU(idx_x,idx_y,:) = I_in;

for ii = 1 : size(I_in,3)
    t = I1LRU(:,:,ii);
    t = imfilter(t',BaseCoeff_x,flag_edge);
    I1LRU(:,:,ii) = imfilter(t',BaseCoeff_y,flag_edge); 
end

I_in=I1LRU;

end