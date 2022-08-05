function [ I_out ] = interpCascade_offcenter( I_in,ratio,edge,window,flag_edge )
%INTERPCASCADE: Vivone's interp23tap using cascade
%   I_in:      input image
%   ratio:     scale ratio between interpolated and input image
%   edge:      edge to remove
%   window:    Window type for FIR filter ('Hamming','Rectangular','Blackman')
%   flag_edge: Matlab boundary options ('circular','symmetric','replicate')

if nargin<=2, edge=0; end
if nargin<=3 || isempty(window), window='Hamming'; end % window='Rectangular';
if nargin<=4 || isempty(flag_edge), flag_edge='circular'; end
if isequal(size(ratio),[2,1]), ratio=ratio.'; end

edge=repmat(edge,3-size(edge));
if numel(ratio)==1, ratio=[ratio,ratio]; end

tol=1e-10;
if any(rem(size(I_in,1)*ratio,1)>tol) || any(rem(size(I_in,2)*ratio,1)>tol)
    error('Target image''s dimensions are not integer');
end

ratio_den=[1,1];
while rem(ratio(1)*ratio_den(1),1)>tol && ratio_den(1)<20
    ratio_den(1)=ratio_den(1)+1;
end
while rem(ratio(2)*ratio_den(2),1)>tol && ratio_den(2)<20
    ratio_den(2)=ratio_den(2)+1;
end

if any(ratio_den>=20)
    error('Ratio is not fractional');
end
ratio_num=round(ratio.*ratio_den);
fact_ratio1=fliplr(factor(ratio_num(1)));
fact_ratio2=fliplr(factor(ratio_num(2)));
fact_ratio1=padarray(fact_ratio1,[0,max(0,length(fact_ratio2)-length(fact_ratio1))],1,'post');
fact_ratio2=padarray(fact_ratio2,[0,max(0,length(fact_ratio1)-length(fact_ratio2))],1,'post');
fact_ratio=[fact_ratio1;fact_ratio2]';

first=[1,1];
for ii=1:size(fact_ratio,1)
    curr_ratio=fact_ratio(ii,:);
    tap=11*curr_ratio+rem(curr_ratio,2);
    edge_extension=repmat(tap,[2,1])/2-edge;
    I_in=edge_extend(I_in,edge_extension);
    edge=(edge+edge_extension).*repmat(curr_ratio,[2,1]);
    I_in=interp23tapGeneral_offcenter(I_in,curr_ratio,tap,window,flag_edge,first);
    if curr_ratio(1)==2, first(1)=0; end
    if curr_ratio(2)==2, first(2)=0; end
    I_out=I_in(edge(1,1)+1:end-edge(2,1),edge(1,2)+1:end-edge(2,2),:);
    if ~isequal(ratio_den,[1,1])
        I_out=imresize(I_out,[size(I_out,1),size(I_out,2)]./ratio_den);
    end
end
