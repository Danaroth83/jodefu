function D = merging_unmixing_MAX(Nodei,Nodej)

global mapwhed

% Specifies the region model of the merging of two regions Ri and Rj

% field size
Si = Nodei.descriptors.size;
Sj = Nodej.descriptors.size;
D.size = Si+Sj;

%field model
leavesi = Nodei.nodeinfo.leaves;
if isempty(leavesi)
    leavesi = Nodei.label;
end
leavesj = Nodej.nodeinfo.leaves;
if isempty(leavesj)
    leavesj = Nodej.label;
end
leaves = union(leavesi,leavesj);
pixels = getpixels(mapwhed,leaves)';
D.size = size(pixels,2);
D.model = unmixing(pixels,'max');

%field boundingbox
bbi = Nodei.descriptors.boundingbox;
bbj = Nodej.descriptors.boundingbox;
D.boundingbox = [min(bbi(1),bbj(1)) max(bbi(2),bbj(2)) ...
    min(bbi(3),bbj(3)) max(bbi(4),bbj(4))];