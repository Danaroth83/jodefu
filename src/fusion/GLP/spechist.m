%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function name : spechist
%
% Inputs :	im = image whos histogram is to be modified
%		h1 = Reference Histogram
%       mbn = maximum bin numer
%
% Outputs :	fin = Histogram Matched Image
%
% Created by : Mr. Dinu Coltuc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fin=spechist(im,h1,mbn)


f1=[0 1 0;
    1 0 1;
	 0 1 0]/4;

f2=[1 0 1;
    0 0 0;
    1 0 1]/4;
 
f3=[0 0 1 0 0;
    0 0 0 0 0;
    1 0 0 0 1
    0 0 0 0 0
 	 0 0 1 0 0]/4;
 
f4=[0 1 0 1 0;
    1 0 0 0 1;
    0 0 0 0 0
    1 0 0 0 1
    0 1 0 1 0]/8;
 
f5=[1 0 0 0 1;
    0 0 0 0 0;
    0 0 0 0 0
    0 0 0 0 0
 	 1 0 0 0 1]/4;
 
[N,M]=size(im);
mas=ones(N,M);
h2=h1;

%L=256;%The only change that I made to make this function work for 16 bit data
%L=65536;
L=mbn;

imag=im*L^5+conv2(im,f1,'same')*L^4+conv2(im,f2,'same')*L^3+conv2(im,f3,'same')*L^2+conv2(im,f4,'same')*L+conv2(im,f5,'same');
sir=reshape(imag,1,M*N);
[ord,ind]=sort(sir);


fin=zeros(N,M);

p1=1; 
p2=1;

for pix=1:N*M,
      x=ind(pix);
      i=rem(x-1,N);
      j=(x-1-i)/N;
      if mas(i+1,j+1)>0
         if h2(p2)==0
            while h2(p2)==0,
               p2=p2+1;
            end
         end         
         fin(i+1,j+1)=p2-1;
         h2(p2)=h2(p2)-1;
      else
         if h1(p1)==0
            while h1(p1)==0,
               p1=p1+1;
            end
         end
         fin(i+1,j+1)=p1-1;
         h1(p1)=h1(p1)-1;
      end         
end
   

