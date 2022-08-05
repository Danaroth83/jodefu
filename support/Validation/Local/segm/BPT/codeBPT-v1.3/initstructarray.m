function tree = initstructarray(handleR,handleO)

fprintf('Initialization of the tree\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

global mapwhed

[m n] = size(mapwhed);
N = max(mapwhed(:));
tree(2*N-1) = struct;
boxes = regionprops(mapwhed,'Boundingbox');
% keyboard
for i=1:N
    tree(i).label = i;
    tree(i).descriptors = handleR(i);
    bb = getboundingbox(boxes(i),[m n]);
    tree(i).descriptors.boundingbox = bb;
    tree(i).nodeinfo = initnodeinfo;
    tree(i).construction.neighbors = getneighlabel(mapwhed(bb(1):bb(2),bb(3):bb(4)),i);
    tree(i).construction.alldist = zeros(size(tree(i).construction.neighbors));
end
fprintf('Creating leaf structures: done\n')
% keyboard
for i=1:N
    neigh = tree(i).construction.neighbors;
    N = neigh(neigh>tree(i).label);
    D = zeros(size(N));
    for j=1:length(N)
        D(j) = handleO(tree(i).descriptors.model,tree(N(j)).descriptors.model);
        ind = tree(N(j)).construction.neighbors == i;
        tree(N(j)).construction.alldist(ind) = D(j);
    end
    tree(i).construction.alldist(neigh>tree(i).label) = D;
    tree(i).construction.dist = min(tree(i).construction.alldist);
    ind = find(tree(i).construction.alldist==min(tree(i).construction.alldist),1,'first');
    tree(i).nodeinfo.sibling = neigh(ind);
    tree(i).nodeinfo.leaves = [];
    tree(i).nodeinfo.iteration = 0;
end
fprintf('Updating leaf structures: done\n')
