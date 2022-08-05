function X = getpixels(map,lbl)

global I

Xc = find(ismember(map,lbl));
d = size(I,2);
X = zeros(length(Xc),d);
for i=1:length(Xc)
    X(i,:) = I(Xc(i),:);
end