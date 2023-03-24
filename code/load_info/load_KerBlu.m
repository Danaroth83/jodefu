function [ KerBlu, start_pos ] = load_KerBlu( im_tag, ratio, GNyq_MS, maxlength, banduniform )

%LOAD_KERBLU Loads Blurring Filter and alignment

if nargin<=3 || isempty(maxlength), maxlength=41; end
if nargin<=4, banduniform=true; end
if banduniform==true, GNyq_MS=mean(GNyq_MS); end

Nb=length(GNyq_MS);
if ratio==1, KerBlu=shiftdim(ones(Nb,1),-2); start_pos=[1,1]; return; end

if strncmpi(im_tag,'Moffett',7)
    Lfilter=9;
    sigma = (1/(2*(2.7725887)/ratio^2))^0.5;
    KerBlu = fspecial('gaussian',[Lfilter Lfilter],sigma);
    if Nb>1, KerBlu = repmat(KerBlu,[1,1,Nb]); end
else
    Lfilter=41;
    KerBlu=zeros(Lfilter,Lfilter,Nb);
    for ii=1:Nb
        sigma = sqrt(((Lfilter-1)*(1/2/ratio))^2/(-2*log(GNyq_MS(ii))));
        KerBlu_temp = fspecial('gaussian',Lfilter,sigma);
        KerBlu_temp = KerBlu_temp./max(KerBlu_temp(:));
        KerBlu_temp = fwind1(KerBlu_temp,kaiser(Lfilter));
        KerBlu(:,:,ii)=KerBlu_temp;
    end
end
if strncmpi(im_tag,'Moffett',7)
    start_pos=[1,1];
else
    start_pos=[floor(ratio/2)+1, floor(ratio/2)+1];
end

%% Fixing too long Blurring Kernel sizes
while size(KerBlu,1)>maxlength || size(KerBlu,2)>maxlength
    KerBlu=KerBlu(2:end-1,2:end-1,:);
end

end

