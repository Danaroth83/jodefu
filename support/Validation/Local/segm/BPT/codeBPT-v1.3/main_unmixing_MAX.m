clear all
close all
clc

load ('Data/cuprite');
cuprite = cuprite(101:end,161:end,:);

v = [37,19,12];

imrgb = rescale(cuprite(:,:,v(1)),1);
imrgb(:,:,2) = rescale(cuprite(:,:,v(2)),1);
imrgb(:,:,3) = rescale(cuprite(:,:,v(3)),1);

figure, imshow(imrgb);

[m,n,p] = size(cuprite);
global I
I = reshape(cuprite,m*n,p);

mapp = grad(cuprite,'supremum');
global mapwhed
mapwhed = whed(mapp);

tic
global tree
tree = initstructarray(@R_unmixing_vca_MAX,@O_spectral);
toc
 
tic
updatestructarray(@O_spectral,@merging_unmixing_MAX,@priority_size);
toc

tic
completestructarray;
toc

tic
prunedtree = pruneBPTunmixing('max');
toc

save('cuprite_BPT_unmixing_max','tree','prunedtree');