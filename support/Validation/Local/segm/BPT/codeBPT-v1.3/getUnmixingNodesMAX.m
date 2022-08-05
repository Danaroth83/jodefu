function [rmse,idx] = getUnmixingNodesMAX(prunedtree,lbl)

nodes = lbl;
idx = [];
rmse = 0.0;
while ~isempty(nodes)
    l = nodes(end);
    nodes(end) = [];
    node = prunedtree(l);
    if node.pruning
        idx = union(idx,l);
%         if strcmpi(node.descriptors.model.method,'mean')
%             rmse = Inf;
%         else
            rmse = max(rmse,node.descriptors.model.RMSE_max);
%         end
    else
        nodes = union(nodes,node.nodeinfo.children);
    end
end