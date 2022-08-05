function [rmse,idx] = getUnmixingNodesMEAN(prunedtree,lbl)

nodes = lbl;
idx = [];
%rmse = 0.0;
N_idx = 0; RMSE_idx = 0;
while ~isempty(nodes)
    l = nodes(end);
    nodes(end) = [];
    node = prunedtree(l);
    if node.pruning
        idx = union(idx,l);
        N_idx = N_idx + length(prunedtree(l).descriptors.model.RMSE);
        if strcmpi(node.descriptors.model.method,'mean') && node.descriptors.model.RMSE_avg == 0
            RMSE_idx = Inf;
        else
            RMSE_idx = RMSE_idx + sum(prunedtree(l).descriptors.model.RMSE);
        end
    else
        nodes = union(nodes,node.nodeinfo.children);
    end
end

rmse = RMSE_idx/N_idx;