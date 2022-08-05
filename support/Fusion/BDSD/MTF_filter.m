function [ I_out ] = MTF_filter( I_in,MTF_Nyq,scale,N )

if nargin<4
    N=scale*10+1;
    % N=scale*5+1;
end
I_in=double(I_in);

[L1,L2,Nb]=size(I_in);
I_out=zeros(L1,L2,Nb);
for ii=1:Nb
    sigma=sqrt(-((N-1)/scale/2)^2/2/log(MTF_Nyq(ii)));
    filter=fspecial('gaussian',N,sigma);
    filter=fwind1(filter./max(filter(:)),kaiser(N));
    I_out(:,:,ii)=imfilter(I_in(:,:,ii),real(filter),'replicate');
    % I_out(:,:,ii)=imfilter(I_in(:,:,ii),real(filter),'symmetric');
end
%{
close all;
Nfft=N*10+1;hi=abs(fftshift(fft(fft(filter,Nfft,1),Nfft,2)));n1=-0.5+1/2/Nfft:1/Nfft:0.5-1/2/Nfft;figure(1);plot(n1,hi(Nfft/2+0.5,:)/max(hi(:)))
hold on; plot([1/scale/2,1/scale/2],[0,1]); plot([-0.5,0.5],[MTF_Nyq,MTF_Nyq]);
figure(2); surf(-N/2+0.5:N/2-0.5,-N/2+0.5:N/2-0.5,filter);
%}
end

