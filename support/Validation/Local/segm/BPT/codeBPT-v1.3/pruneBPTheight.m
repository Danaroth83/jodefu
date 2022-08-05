function prunedtree = pruneBPTheight(h)

global tree

prunedtree = tree;
N = length(prunedtree);
Nbleaves = (N+1)/2;
for i=1:N
    prunedtree(i).pruning = 0;
end

H = [prunedtree.nodeinfo];
H = [H.height];

ind = false(1,N);
ind(1:Nbleaves) = H(1:Nbleaves)<=h;
ind(Nbleaves+1:end) = H(Nbleaves+1:end)==h;

node = find(ind);
for i=1:length(node)
    prunedtree(node(i)).pruning = 1;
end