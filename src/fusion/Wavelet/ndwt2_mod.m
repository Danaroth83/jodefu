function ndwt2_mod(I_in, m, filter_b, filter_p)

if nargin<=3 
    filter_p=[-0.0000248,0.003075,-0.04159,0.196,-0.4584,0.6018];
    filter_p=[filter_p,filter_p(end-1:-1:1)];
end
if nargin<=2
    filter_b=[1,4,1]/6;
end

s1=imfilter(I_in,filter_b);
keyboard
end