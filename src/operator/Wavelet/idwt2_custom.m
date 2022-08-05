%% 2D Inverse Discrete Wavelet Transform
%
% Description:
% Performs the Inverse 2D Discrete Wavelet Transform over the first two 
% dimensions of a given image, assuming the products were obtained with
% the dwt2_custom script
%
% Usage:
% y=idwt2_custom(x,lo,hi,lev,sizes)
%
% Input:
% y    : Nd matrix
% lo   : Low Pass Reconstruction Filter
% hi   : High Pass Reconstruction Filter
% lev  : Amount of Decomposition Levels
% sizes: Sizes of the final image (default: Sizes of the input)
% 
% Output:
% y   : Inverse 2D DWT over the first two dimensions

function x=idwt2_custom(x,lo,hi,lev,sizes)

    if nargin<=3, lev=1; end
    if nargin<=4, sizes=size(x); end
    
    sizex=size(x);
    if sizex(2)==1, idwt_custom(x,lo,hi,lev,sizes); return; end
    x=reshape(x,size(x,1),size(x,2),[]);

    L1=size(x,1)/2^(lev-1);
    L2=size(x,2)/2^(lev-1);
    
    for jj=1:lev
        for ii=1:size(x,3)
            ca=x(1:L1/2,1:L2/2,ii);
            ch=x(L1/2+1:L1,1:L2/2,ii);
            cv=x(1:L1/2,L2/2+1:L2,ii);
            cd=x(L1/2+1:L1,L2/2+1:L2,ii);
            x(1:L1,1:L2,ii)=idwt2(ca,ch,cv,cd,lo,hi,'mode','per');
        end
        L1=L1*2;
        L2=L2*2;
    end
    
    modif=2^lev;
    L1=mod(sizes(1),modif);
    L2=mod(sizes(2),modif);
    x=x(floor(L1/2)+1:end-ceil(L1/2),floor(L2/2)+1:end-ceil(L2/2),:);
    
    sizex(1)=size(x,1);
    sizex(2)=size(x,2);
    x=reshape(x,sizex);
    
end