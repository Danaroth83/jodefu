function [ I_out ] = interpCascadeAlparone_Matlabedge( I_in,ratio,edge,flag_edge )
%INTERPCASCADEALPARONE Generalized version of Alparone's interp23tap
%   I_in:      input image
%   ratio:     scale ratio between interpolated and input image
%   edge:      edge to remove
%   flag_edge: Matlab boundary options ('circular','symmetric','replicate')

if nargin<=3 || isempty(flag_edge)
    flag_edge='circular';
end
if nargin<=2 || isempty(edge)
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
        order=11;
        edge=edge*curr_ratio;
        if first==1
            I_in=interp23tap_mod2(I_in,curr_ratio,order,0,flag_edge);
            first=0;
        else
            I_in=interp23tap_mod2(I_in,curr_ratio,order,1,flag_edge);
        end
    else
        order=11;
        edge=edge*curr_ratio;
        I_in=interp23tap_mod2(I_in,curr_ratio,order,0,flag_edge);
    end
    I_out=I_in(edge+1:end-edge,edge+1:end-edge,:);
    if ratio_den~=1
        I_out=imresize(I_out,1/ratio_den);
    end
end

