function msexp = expand4(ms)
%
% msexp = expand4(ms);
%
% keyboard
load h67
Nb = size(ms,3);
msexp = imresize(ms,4);
for n=1:Nb
    msexp(:,:,n) = filter_image(msexp(:,:,n),h67); 
end
