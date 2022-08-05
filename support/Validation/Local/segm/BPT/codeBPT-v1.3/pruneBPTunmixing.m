function prunedtree = pruneBPTunmixing(method)

global tree;

% Initialize
prunedtree = tree;
N = length(prunedtree);
Nbleaves = (N+1)/2;
for i=1:Nbleaves
    prunedtree(i).pruning = true;
    prunedtree(i).isleaf = true;
end
for i=Nbleaves+1:N
    prunedtree(i).pruning = false;
    prunedtree(i).isleaf = false;
end

% create model if need
global mapwhed
if ~isfield(tree(N).descriptors.model,'E')
    for i=1:N
        if i<=Nbleaves
            lbl = i;
        else
            lbl = prunedtree(i).nodeeinfo.leaves;
        end
        pixels = getpixels(mapwhed,lbl)';
        prunedtree(i).descriptors.model = unmixing(pixels);
    end
end

% Prune
if strcmpi(method,'max')
    prunedtree = prunerecursiveBPTunmixingMAX(prunedtree,N);
else
    prunedtree = prunerecursiveBPTunmixingMEAN(prunedtree,N);
end
