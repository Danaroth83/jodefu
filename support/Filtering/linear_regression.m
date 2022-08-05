function [ w ] = linear_regression( I_des, I_in )
%LINEAR_REGRESSION Summary of this function goes here
%   Detailed explanation goes here

[L1,L2,Nb]=size(I_in);
I_des=reshape(I_des,[L1*L2,size(I_des,3)]);
I_in=reshape(I_in,[L1*L2,Nb]);

Rx=I_in'*I_in;
p=I_in'*I_des;
w=Rx\p;

end

