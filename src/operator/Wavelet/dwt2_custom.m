%% 2D Discrete Wavelet Transform
%
% Description:
% Performs the 2D Discrete Wavelet Transform over the first two dimension
% of a given image; the results are spatially concatenated as usually done
% in wavelet literature (LL,HL,LH,HH composition); the algorithm adds edges
% to reach the proper power of 2 divisibility
% Note: If the second dimension is a singleton, 1D DWT is performed instead
%
% Usage:
% y=dwt2_custom(x,lo,hi,lev)
%
% Input:
% x   : Nd matrix
% lo  : Low Pass Decomposition Filter
% hi  : High Pass Decomposition Filter
% lev : Amount of Decomposition Levels
% 
% Output:
% y   : Nd matrix where the DWT was applied over the first 2 dimensions,
%       The edges are extended to reach a multiple of 2^lev

function y=dwt2_custom(x,lo,hi,lev)

    sizex=size(x);
    if sizex(2)==1, y=dwt_custom(x,lo,hi,lev); return; end
    x=reshape(x,size(x,1),size(x,2),[]);

    modif=2^lev;
    L1=mod(-size(x,1),modif);
    L2=mod(-size(x,2),modif);
    x=padarray(x,[ceil(L1/2),ceil(L2/2)],'post');
    x=padarray(x,[floor(L1/2),floor(L2/2)],'pre');

    y=zeros(size(x));
    ca=x;
    for jj=1:lev
        [L1,L2,~]=size(ca);
        ca_new=zeros(L1/2,L2/2,size(x,3));
        for ii=1:size(x,3)
            [ca_new(:,:,ii),ch,cv,cd]=dwt2(ca(:,:,ii),lo,hi,'mode','per');
            y(L1/2+1:L1,1:L2/2,ii)=ch;
            y(1:L1/2,L2/2+1:L2,ii)=cv;
            y(L1/2+1:L1,L2/2+1:L2,ii)=cd;
        end
        ca=ca_new;
    end
    y(1:L1/2,1:L2/2,:)=ca;
        
    sizex(1)=size(y,1);
    sizex(2)=size(y,2);
    y=reshape(y,sizex);
    
end