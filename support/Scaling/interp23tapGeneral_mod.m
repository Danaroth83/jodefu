function I_Interpolated = interp23tapGeneral_mod(I_Interpolated,ratio,L,flag_align,window,flag_edge)
% keyboard
if nargin<=5 || isempty(flag_edge)
    flag_edge='circular';  % Also available 'symmetric'
end
if nargin<=4 || isempty(window)
    window='Hamming';
end
if nargin<=3 || isempty(flag_align)
    flag_align=0;
end
if nargin<=2 || isempty(L)
    L = 44; % tap
end

[r,c,b] = size(I_Interpolated);

switch window
    case 'Rectangular'
        BaseCoeff = ratio.*fir1(L,1./ratio,rectwin(L+1),'noscale');
    case 'Blackman'
        BaseCoeff = ratio.*fir1(L,1./ratio,blackman(L+1),'noscale');
    otherwise
        BaseCoeff = ratio.*fir1(L,1./ratio,'noscale');
end

I1LRU = zeros(ratio.*r, ratio.*c, b); 
if flag_align==0 || rem(ratio,2)~=0
    I1LRU(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:,:) = I_Interpolated;
else
    I1LRU(floor(ratio/2):ratio:end,floor(ratio/2):ratio:end,:,:) = I_Interpolated;
end

for ii = 1 : b
    t = I1LRU(:,:,ii);
    t = imfilter(t',BaseCoeff,flag_edge); 
    I1LRU(:,:,ii) = imfilter(t',BaseCoeff,flag_edge); 
end

I_Interpolated = I1LRU;

end