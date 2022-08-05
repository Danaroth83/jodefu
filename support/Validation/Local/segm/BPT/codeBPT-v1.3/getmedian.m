function med = getmedian(map)

global I

med = zeros(max(map(:)),size(I,2));
for i=1:max(map(:))
    pix = getpixels(map,i);
    med(i,:) = median(pix);
end