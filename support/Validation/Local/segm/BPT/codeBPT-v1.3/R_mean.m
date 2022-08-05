function D = R_mean(lbl)

% Creates a structure. The region model of the set of pixels pixls is its
% mean spectrum

global mapwhed
pixels = getpixels(mapwhed,lbl)';

D = struct;
D.size = size(pixels,2);
D.model = mean(pixels,2);