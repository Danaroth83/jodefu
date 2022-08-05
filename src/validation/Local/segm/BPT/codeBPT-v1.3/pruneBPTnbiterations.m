function prunedtree = pruneBPTnbiterations(NbIter)

global tree
% keyboard
prunedtree = tree;
N = length(prunedtree);
Nbleaves = (N+1)/2;
for i=1:N
    prunedtree(i).pruning = 0;
end

IterMax = prunedtree(end).nodeinfo.iteration;
NbIter = round(NbIter);

if NbIter>=IterMax
    prunedtree(end).pruning = 1;
else
    IterCut = NbIter;
    SetLeaves = 1:Nbleaves;
    
    SetNode = [];
    while ~isempty(SetLeaves)
        leaf = SetLeaves(1);
        p = prunedtree(leaf).nodeinfo.branch;
        IterBranch = [prunedtree(p).nodeinfo];
        IterBranch = [IterBranch.iteration];
        node = p(find(IterBranch<=IterCut,1,'last'));
        if isempty(node)
            node = leaf;
            SetLeaves(1) = [];
        else
            leaves = prunedtree(node).nodeinfo.leaves;
            SetLeaves = setdiff(SetLeaves,leaves);
        end
        SetNode = cat(2,SetNode,node);
    end
    
    for i=1:length(SetNode)
        prunedtree(SetNode(i)).pruning = 1;
    end
    
end

