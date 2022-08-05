%%%_______________________________________________________________
%%%
function out = filter_image(in,filter)
%
% out = filter_image(in,filter);
%
% keyboard
[dimy dimx Nb] = size(in);

h = filter;
H = h'*h;
[winy winx] = size(H);

wx = (winx-1)/2;
wy = (winy-1)/2;
ext = zeros(dimy+winy-1,dimx+winx-1);
ext(wy+1:wy+dimy,wx+1:wx+dimx) = in;

for k = 1:wy
   ext(wy-k+1,wx+1:wx+dimx) = in(k+1,:);
   ext(wy+dimy+k,wx+1:wx+dimx) = in(dimy-k,:);
end
for k = 1:wx
   ext(:,wx-k+1) = ext(:,wx+k+1);
   ext(:,wx+dimx+k) = ext(:,dimx+wx-k);
end

out = zeros(dimy,dimx,Nb);
for i=1:Nb
    filt_full = imfilter(ext(:,:,i),H,'full');
    out(:,:,i) = filt_full(winy:winy+dimy-1,winx:winx+dimx-1);
end
