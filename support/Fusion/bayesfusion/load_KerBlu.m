function [ KerBlu, start_pos ] = load_KerBlu( im_tag, ratio, GNyq_MS, maxlength )
%LOAD_KERBLU Loads Blurring Filter and alignment
if nargin<=3
    maxlength=41;
end

if strncmpi(im_tag,'Moffett',7)
    Lfilter=9;
    sigma = (1/(2*(2.7725887)/ratio^2))^0.5;
    KerBlu = fspecial('gaussian',[Lfilter Lfilter],sigma);
else
    Lfilter=41;
    sigma = sqrt(((Lfilter-1)*(1/2/ratio))^2/(-2*log(mean(GNyq_MS))));
    KerBlu = fspecial('gaussian',Lfilter,sigma);
    KerBlu = KerBlu./max(KerBlu(:));
    KerBlu = fwind1(KerBlu,kaiser(Lfilter));
end
if strncmpi(im_tag,'Moffett',7)
    start_pos=[1,1];
else
    start_pos=[floor(ratio/2)+1, floor(ratio/2)+1];
end

%% Fixing too long Blurring Kernel sizes
while size(KerBlu,1)>maxlength || size(KerBlu,2)>maxlength
    KerBlu=KerBlu(2:end-1,2:end-1);
end

end

