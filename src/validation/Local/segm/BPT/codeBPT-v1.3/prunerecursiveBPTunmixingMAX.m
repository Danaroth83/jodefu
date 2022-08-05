function prunedtree = prunerecursiveBPTunmixingMAX(prunedtree,root_lbl)

% subtrees
children = prunedtree(root_lbl).nodeinfo.children;
child1 = prunedtree(children(1));
child2 = prunedtree(children(2));

% get the indices and the max RMSE values of the subtrees
if child1.isleaf
    if strcmpi(child1.descriptors.model.method,'mean') && child1.descriptors.model.RMSE_max == 0
        rmse1 = Inf;
    else
        rmse1 = child1.descriptors.model.RMSE_max;
    end
    idx1 = child1.label;
else
    prunedtree = prunerecursiveBPTunmixingMAX(prunedtree,child1.label);
    [rmse1,idx1] = getUnmixingNodesMAX(prunedtree,child1.label);
end
if child2.isleaf
    if strcmpi(child2.descriptors.model.method,'mean') && child2.descriptors.model.RMSE_max == 0
        rmse2 = Inf;
    else
        rmse2 = child2.descriptors.model.RMSE_max;
    end
    idx2 = child2.label;
else
    prunedtree = prunerecursiveBPTunmixingMAX(prunedtree,child2.label);
    [rmse2,idx2] = getUnmixingNodesMAX(prunedtree,child2.label);
end

% Calculate max RMSE
if prunedtree(root_lbl).descriptors.model.RMSE_max <= max(rmse1,rmse2)
    prunedtree(root_lbl).pruning = true;
    for i=1:length(idx1)
        prunedtree(idx1(i)).pruning = false;
    end
    for i=1:length(idx2)
        prunedtree(idx2(i)).pruning = false;
    end
end