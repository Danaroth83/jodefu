function I_Interpolated = interp23tapGeneral_mod2(I_Interpolated,ratio,L,flag_align,window,flag_edge)
% keyboard
if nargin<=5 || isempty(flag_edge), flag_edge='circular'; end % Also available 'symmetric'
if nargin<=4 || isempty(window), window='Hamming'; end
if nargin<=3 || isempty(flag_align), flag_align=0; end
if nargin<=2 || isempty(L), L = 44; end % tap
if numel(ratio==1), ratio=[ratio,ratio]; end

[r,c,b] = size(I_Interpolated);

switch window
    case 'Rectangular'
        BaseCoeff_x = ratio(1).*fir1(L,1./ratio(1),rectwin(L+1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L,1./ratio(2),rectwin(L+1),'noscale');
    case 'Blackman'
        BaseCoeff_x = ratio(1).*fir1(L,1./ratio(1),blackman(L+1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L,1./ratio(2),blackman(L+1),'noscale');
    otherwise
        BaseCoeff_x = ratio(1).*fir1(L,1./ratio(1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L,1./ratio(2),'noscale');
end

I1LRU = zeros(ratio(1).*r, ratio(2).*c, b);

if flag_align==0 || rem(ratio(1),2)~=0
    idx_x=floor(ratio(1)/2)+1:ratio(1):size(I1LRU,1);
else
    idx_x=floor(ratio(1)/2):ratio(1):size(I1LRU,1);
end
if flag_align==0 || rem(ratio(1),2)~=0
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