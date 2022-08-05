function D = R_unmixing_vca_MEAN(lbl)

% Creates a structure. The region model of the set of pixels is given by
% the set of induced endmembers by VCA method

global mapwhed
pixels = getpixels(mapwhed,lbl)';
D = struct;
D.size = size(pixels,2);
D.model = unmixing(pixels,'mean');
