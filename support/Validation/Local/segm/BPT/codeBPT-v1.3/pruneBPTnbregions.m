function prunedtree = pruneBPTnbregions(NbReg)
% keyboard
global tree

prunedtree = tree;
N = length(prunedtree);
Nbleaves = (N+1)/2;
for i=1:N
    prunedtree(i).pruning = 0;
end

NbReg = round(NbReg);

if NbReg>Nbleaves
    
    for i=1:Nbleaves
        prunedtree(i).pruning = 1;
    end
    
elseif NbReg>0
    
    IterMax = prunedtree(end).nodeinfo.iteration;
    prunedtree = pruneBPTnbiterations(IterMax-NbReg+1);
    
else
    error('desired number of regions must be a positive integer')    
end