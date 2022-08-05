function I_in = interp23tapGeneral_offcenter(I_in,ratio,L,window,flag_edge,flag_alignright)
% keyboard
if nargin<=2 || isempty(L), L = 44; end % tap
if nargin<=3 || isempty(window), window='Hamming'; end
if nargin<=4 || isempty(flag_edge), flag_edge='circular'; end % Also available 'symmetric'
if nargin<=5 || isempty(flag_alignright), flag_alignright=1; end
if numel(ratio)==1, ratio=[ratio,ratio]; end
if numel(L)==1, L=[L,L]; end
if numel(flag_alignright)==1, flag_alignright=[flag_alignright,flag_alignright]; end

switch window
    case 'Rectangular'
        BaseCoeff_x = ratio(1).*fir1(L(1),1./ratio(1),rectwin(L(1)+1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L(2),1./ratio(2),rectwin(L(2)+1),'noscale');
    case 'Blackman'
        BaseCoeff_x = ratio(1).*fir1(L(1),1./ratio(1),blackman(L(1)+1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L(2),1./ratio(2),blackman(L(2)+1),'noscale');
    otherwise
        BaseCoeff_x = ratio(1).*fir1(L(1),1./ratio(1),'noscale');
        BaseCoeff_y = ratio(2).*fir1(L(2),1./ratio(2),'noscale');
end

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

for ii = 1 : b
    t = I1LRU(:,:,ii);
    t = imfilter(t',BaseCoeff_x,flag_edge); 
    I1LRU(:,:,ii) = imfilter(t',BaseCoeff_y,flag_edge); 
end

I_in=I1LRU;

end