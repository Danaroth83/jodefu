function [ I_out ] = interpCascade_customsymmetric( I_in,ratio,edge,window )
%INTERPCASCADE Generalized version of interp23tap using cascade
% This is a custom version of MATlab 'symmetric' imfilter boundary option
%   I_in:    input image
%   ratio:   scale ratio between interpolated and input image
%   edge:    edge to remove
%   window:  Window type for FIR filter ('Hamming','Rectangular','Blackman')
if nargin<=3
    window='Hamming';
end
if nargin<=2
    edge=0;
end

tol=0.0001;
if mod(size(I_in,1)*ratio,1)>tol || mod(size(I_in,1)*ratio,1)>tol
    error('Target image''s dimensions are not integer');
end

ratio_den=1;
while mod(ratio*ratio_den,1)>tol && ratio_den<20
    ratio_den=ratio_den+1;
end
if ratio_den>=20
    error('Ratio is not fractional');
end
ratio_num=ratio*ratio_den;
fact_ratio=fliplr(factor(ratio_num));

first=1;
for ii=1:length(fact_ratio)
    curr_ratio=fact_ratio(ii);
    if rem(curr_ratio,2)==0
        tap=11*curr_ratio;
        edge_extension=tap/2-edge;
        if edge_extension>0
            I_in=edge_extend(I_in,edge_extension,edge_extension,'r');
        elseif edge_extension<0
            I_in=I_in(-edge_extension+1:end+edge_extension,-edge_extension+1:end+edge_extension,:);
        end
        edge=(edge+edge_extension)*curr_ratio;
        if first==1
            I_in=interp23tapGeneral_mod(I_in,curr_ratio,tap,0,window);
            first=0;
        else
            I_in=interp23tapGeneral_mod(I_in,curr_ratio,tap,1,window);
        end
    else
        tap=11*curr_ratio+1;
        edge_extension=tap/2-edge;
        if edge_extension>0
            I_in=edge_extend(I_in,edge_extension,edge_extension,'r');
        elseif edge_extension<0
            I_in=I_in(-edge_extension+1:end+edge_extension,-edge_extension+1:end+edge_extension,:);
        end
        edge=(edge+edge_extension)*curr_ratio;
        I_in=interp23tapGeneral_mod(I_in,curr_ratio,tap,0,window);
    end
    I_out=I_in(edge+1:end-edge,edge+1:end-edge,:);
    if ratio_den~=1
        I_out=imresize(I_out,1/ratio_den);
    end
end

