function prunedtree = prunerecursiveBPTunmixingMEAN(prunedtree,root_lbl)

% subtrees
children = prunedtree(root_lbl).nodeinfo.children;
child1 = prunedtree(children(1));
child2 = prunedtree(children(2));

% get the indices of the subtrees
if child1.isleaf
    idx1 = child1.label;
else
    prunedtree = prunerecursiveBPTunmixingMEAN(prunedtree,child1.label);
    [~,idx1] = getUnmixingNodesMEAN(prunedtree,child1.label);
end
if child2.isleaf
    idx2 = child2.label;
else
    prunedtree = prunerecursiveBPTunmixingMEAN(prunedtree,child2.label);
    [~,idx2] = getUnmixingNodesMEAN(prunedtree,child2.label);
end
idx = union(idx1,idx2);

% calculate average RMSE
N_idx = 0; RMSE_idx = 0;
for i=1:length(idx)
    N_idx = N_idx + length(prunedtree(idx(i)).descriptors.model.RMSE);
    if strcmpi(prunedtree(idx(i)).descriptors.model.method,'mean') && prunedtree(idx(i)).descriptors.model.RMSE_avg == 0
        RMSE_idx = Inf;
    else
        RMSE_idx = RMSE_idx + sum(prunedtree(idx(i)).descriptors.model.RMSE);
    end
end
RMSE_idx = RMSE_idx/N_idx;
if prunedtree(root_lbl).descriptors.model.RMSE_avg <= RMSE_idx;
    prunedtree(root_lbl).pruning = true;
    for i=1:length(idx1)
        prunedtree(idx1(i)).pruning = false;
    end
    for i=1:length(idx2)
        prunedtree(idx2(i)).pruning = false;
    end
end