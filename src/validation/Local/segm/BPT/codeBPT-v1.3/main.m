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

[m n p] = size(cuprite);
global I
I = reshape(cuprite,m*n,p);

global mapp
mapp = grad(cuprite,'supremum');
global mapwhed
mapwhed = whed(mapp);

%% creation of the tree
tic
global tree
tree = initstructarray(@R_mean,@O_SAM);
toc

tic
updatestructarray(@O_SAM,@merging_mean,@priority_size);
toc

completestructarray;

%% pruning of the tree

% prunedtree = pruneBPTheight(40);
% prunedtree = pruneBPTnbiterations(2850);
prunedtree = pruneBPTnbregions(25);
retrievesegmentation(prunedtree,imrgb);
