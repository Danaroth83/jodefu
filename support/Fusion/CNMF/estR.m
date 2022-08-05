function R = estR(HS,MS,mask)
%--------------------------------------------------------------------------
% Estimation of relative SRFs
%
% USAGE
%       R = estR(HS,MS,mask)
%
% INPUT
%       HS  : Low-spatial-resolution HS image (rows2,cols2,bands2)
%       MS  : MS image (rows1,cols1,bands1)
%       mask: (optional) Binary mask for processing (rows2,cols2)
%
% OUTPUT
%       R   : Relative SRFs (bands1,bands2)
%
% Copyright (c) 2015 Naoto Yokoya All rights reserved.
% Email: yokoya@sal.rcast.u-tokyo.ac.jp
% Update: 2015/02/10
%
% REFERENCE
% N. Yokoya, N. Mayumi, and A. Iwasaki, "Cross-calibration for data 
% fusion of EO-1/Hyperion and Terra/ASTER," IEEE J. Sel. Topics Appl. 
% Earth Observ. Remote Sens., vol. 6, no. 2, pp. 419-426, 2013.
%--------------------------------------------------------------------------

if nargin == 2
    masking = 0;
elseif nargin == 3
    masking = 1;
end

[rows1,cols1,bands1] = size(MS);
[rows2,cols2,bands2] = size(HS);
if masking == 1
    HS = reshape([reshape(HS,[],bands2) reshape(mask,[],1)],rows2,cols2,[]);
    bands2 = size(HS,3);
end
R = zeros(bands1,bands2);

% downgrade spatial resolution
LR_MS = zeros(rows2,cols2,bands1);
for b = 1:bands1
    LR_MS(:,:,b) = imresize(reshape(MS(:,:,b),rows1,cols1),[rows2 cols2]);
end
A = reshape(HS,[],bands2);
H = A'*A;

options = optimset('Display','off','MaxIter',500);

for k = 1:bands1
    b = double(reshape(LR_MS(:,:,k),[],1));
    f = -A'*b;
    C = -eye(bands2);
    e = zeros(bands2,1);    
    x = quadprog(H,f,C,e,[],[],[],[],[],options);
    R(k,:) = reshape(x,1,[]);
end
